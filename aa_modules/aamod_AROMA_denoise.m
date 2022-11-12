function [ aap,resp ] = aamod_AROMA_denoise(aap, task, subject_index, session_index)
%
% denoise an EPI using FSL's ICA-AROMA
%
% This requires AROMA-ICA to be installed and added as an aa toolbox. The commit of e8d7a5891563aff7ff71375078314c417d32d808 is recommended.
% Also, AROMA requires Python 3.6 and pip.
%
% The easiest way to install is prolly:
%
%   % cd /users/abcd1234/tools
%   % git clone https://github.com/maartenmennes/ICA-AROMA.git
%   % cd ICA-AROMA
%   % git reset --hard e8d7a5891563aff7ff71375078314c417d32d808
%   % pip install -r requirements.txt
%
%   (assuming that the repo link is still valid)
%
% add the correspodning entry to your parameterset
%   <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
%       <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>aroma</name>
%           <dir ui='dir'>/users/abcd1234/tools/ICA-AROMA</dir>
%   </toolbox>
%
% note this is unsupervised denoising, which as a rule does not
% perform as well as using ICA with training data...
%
% Change History
%
% summer 2021 [MSJ] - new
%

resp='';

switch task
	
    case 'report'
				
    case 'doit'
        
        working_dir = aas_getsesspath(aap, subject_index, session_index);
        
        save_fsloutputtype = aap.directory_conventions.fsloutputtype;
        aap.directory_conventions.fsloutputtype = 'NIFTI_GZ';
        % NB: AROMA seems to require full paths to everything
        [~, AROMA] = aas_cache_get(aap,'aroma');
        AROMA_FNAME = fullfile(AROMA.toolPath,'ICA_AROMA.py');
     
        % we assume EPI and structural are in native space
        % we assume epi is realigned and coregistered to structural
        % (although we will acually warm epi to FSL MNI template and run
        % AROMA in standard space)
        
        EPI_FNAME = aas_getimages_bystream(aap, subject_index, session_index, 'epi');
        STRUCTURAL_FNAME = aas_getfiles_bystream(aap, subject_index, 'structural');

        % Previously we used mcflirtto generate .par input for AROMA
        % Maartin (FSL) suggested using SPM rp.txt file generated by SPM realign 
        % (it was noted SPM rp.txt and FSL .par files have different ordering but that's okay. Go figure).
         
        MC_FNAME = aas_getimages_bystream(aap, subject_index, session_index, 'realignment_parameter');
          
        % this is the MNI structural reference. FSL comes with a number of
        % them preinstalled. We use the 2mm, because fnirt takes > 5
        % hours to run using the 1mm and 30 min using 2mm.
             
        REFERENCE_FNAME = fullfile(aap.directory_conventions.fsldir,'data/standard/MNI152_T1_2mm_brain.nii.gz');
        REFMASK_FNAME = fullfile(aap.directory_conventions.fsldir,'data/standard/MNI152_T1_2mm_brain_mask.nii.gz');

        % output will go in WD/DENOISED/denoised_func_data_nonaggr.nii[.gz]
        % for the vanilla AROMA settings used here
        
        OUT_DIR = fullfile(working_dir,'DENOISE');
      
        % 0) brain extraction
        %
        % bet my_structural my_betted_structural
      
        BETTEDSTRUCTURAL_FNAME = fullfile(working_dir,'betted_structural.nii.gz');
        aas_log(aap, false, sprintf('Brain extraction on %s...', STRUCTURAL_FNAME));      
        command = sprintf('bet %s %s', STRUCTURAL_FNAME, BETTEDSTRUCTURAL_FNAME);
        [ status,result ] = aas_runfslcommand(aap, command);        
        if (status > 0)
           aas_log(aap, true, sprintf('Error running bet: %s', result));
        end     
        
        % 1) generate affmat (structural to MNI)
        %
        % flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in my_betted_structural -omat my_affine_transf.mat
        %
        
        STRUCT2MNI_AFFINE = fullfile(working_dir,'struct2MNI.mat');
        aas_log(aap, false, sprintf('Running FLIRT on %s...', BETTEDSTRUCTURAL_FNAME));      
        command = sprintf('flirt -in %s -ref %s -omat %s', BETTEDSTRUCTURAL_FNAME, REFERENCE_FNAME, STRUCT2MNI_AFFINE);
        [ status,result ] = aas_runfslcommand(aap, command);        
        if (status > 0)
           aas_log(aap, true, sprintf('Error running flirt: %s', result));
        end
        
        % 2) generate warp field (structural to MNI)
        %
        % fnirt --in=my_structural --aff=my_affine_transf.mat --cout=my_nonlinear_transf --config=T1_2_MNI152_2mm
        %
        % note fnirt has a different syntax for specifying options than all the other commands
        %
        
        STRUCT2MNI_WARP = fullfile(working_dir,'struct2MNI_warp.nii.gz');
        aas_log(aap, false, sprintf('Running FNIRT on %s...', STRUCTURAL_FNAME));
        command = sprintf('fnirt --in=%s --aff=%s --cout=%s --config=T1_2_MNI152_2mm', STRUCTURAL_FNAME, STRUCT2MNI_AFFINE, STRUCT2MNI_WARP);
        [ status,result ] = aas_runfslcommand(aap,command);        
         if (status > 0)
            aas_log(aap, true, sprintf('Error running fnirt: %s', result));
         end    
         
        % 3) func-to-struct affine (used in func-to-MNI)
        %
        % flirt -ref my_betted_structural -in my_functional -dof 6 -omat func2struct.mat

        EPI2STRUCT_AFFINE = fullfile(working_dir,'epi2struct.mat');
        aas_log(aap, false, sprintf('Running FLIRT on %s...', EPI_FNAME));
        command = sprintf('flirt -in %s -ref %s -dof 6 -omat %s',EPI_FNAME, BETTEDSTRUCTURAL_FNAME, EPI2STRUCT_AFFINE);
        [ status,result ] = aas_runfslcommand(aap,command);        
        if (status > 0)
            aas_log(aap, true, sprintf('Error running flirt: %s', result));
        end
           
        % generate mean epi (need for detrending)

        MEANEPI_FNAME = fullfile(working_dir, 'meanepi.nii.gz');
        aas_log(aap, false, sprintf('creating mean epi...'));
        command = sprintf('fslmaths %s -Tmean %s', EPI_FNAME, MEANEPI_FNAME);
        [ status,result ] = aas_runfslcommand(aap,command);
        if (status > 0)
           aas_log(aap, true, sprintf('Error creating mean epi: %s', result));
        end

        % DETREND - do a highpass filter based on 100s

        DICOMHEADERS = load(aas_getfiles_bystream(aap,subject_index, session_index,'epi_dicom_header'));
        TR = DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
        bptf = 0.5 * (100/TR);

        DETREND_FNAME = fullfile(working_dir, 'epi_detrend.nii.gz');
        % NB -- need to add meanepi back in after detrend
        aas_log(aap, false, sprintf('detrending %s...', EPI_FNAME));
        command = sprintf('fslmaths %s -bptf %f -1 -add %s %s', EPI_FNAME, bptf, MEANEPI_FNAME, DETREND_FNAME);
        [ status,result ] = aas_runfslcommand(aap,command);
        if (status > 0)
           aas_log(aap, true, sprintf('Error detrending: %s', result));
        end

        SMOOTHED_FNAME = fullfile(working_dir,'epi_smooth.nii.gz');
        aas_log(aap, false, sprintf('smoothing %s...', DETREND_FNAME));
        command = sprintf('fslmaths %s -kernel gauss 2.1233226 -fmean %s', DETREND_FNAME, SMOOTHED_FNAME);
        [ status,result ] = aas_runfslcommand(aap,command);
        if (status > 0)
           aas_log(aap, true, sprintf('Error running smoothing: %s', result));
        end
        
        % 4) func-to-MNI
        %
        % applywarp --in=my_functional --ref=${FSLDIR}/data/standard/MNI152_T1_2mm --warp=my_nonlinear_transf --premat=func2struct.mat --out=my_warped_functional         

        EPIMNI_FNAME = fullfile(working_dir,'epiMNI.nii.gz');        
        WARPREF = fullfile(aap.directory_conventions.fsldir,'data/standard/MNI152_T1_2mm');
        aas_log(aap, false, sprintf('transforming %s to MNI space...', SMOOTHED_FNAME));
        command = sprintf('applywarp --in=%s --ref=%s --warp=%s --premat=%s --out=%s', SMOOTHED_FNAME, WARPREF, STRUCT2MNI_WARP, EPI2STRUCT_AFFINE, EPIMNI_FNAME);
        [ status,result ] = aas_runfslcommand(aap,command);
        if (status > 0)
            aas_log(aap, true, sprintf('Error running applywarp: %s', result));
        end      
                      
        % now finally ready to run AROMA
        
        aas_log(aap, false, sprintf('Running AROMA on %s...', EPIMNI_FNAME));
        
        % AROMA will abort if OUT_DIR already exists (it might from aborted
        % previous run). If we get to here, assume we're cool with deleting it...
        
        if exist(OUT_DIR,'dir')
            rmdir(OUT_DIR,'s');
        end
        
        command = sprintf('python %s -in %s -out %s -mc %s -m %s', AROMA_FNAME, EPIMNI_FNAME, OUT_DIR, MC_FNAME, REFMASK_FNAME);
        [ status,result ] = aas_runfslcommand(aap,command);
        if (status > 0)
            aas_log(aap, true, sprintf('Error running AROMA: %s', result));
        else
            aas_log(aap, false, result);
        end
        
        % desc
        
        % currently, AROMA generates intermediate gzip-ped files, which
        % means the output is also .gz. However, .gz doesn't play nice with
        % SPM and other aa modules, so gunzip here before proceeding.
        % also delete the .gz (bc big)

        output_fname = fullfile(OUT_DIR,'denoised_func_data_nonaggr.nii');
        gunzip([output_fname '.gz'],working_dir);
        delete([output_fname '.gz']);
        aap = aas_desc_outputs(aap, subject_index, session_index, 'epi', spm_file(output_fname,'path',working_dir));
       
        % cleanup 

        if (aap.tasklist.currenttask.settings.deletegzips)
            delete(BETTEDSTRUCTURAL_FNAME);
            delete(STRUCT2MNI_AFFINE);
            delete(STRUCT2MNI_WARP);
            delete(EPI2STRUCT_AFFINE);
            delete(MEANEPI_FNAME);
            delete(DETREND_FNAME);
            delete(SMOOTHED_FNAME);
            delete(EPIMNI_FNAME);
        end
        
        aap.directory_conventions.fsloutputtype = save_fsloutputtype;
        % aside: if you delete DENOISE, it would remove about about 99%
        % of unnecessary files (which might be more than 1 GB!). Maybe make
        % this an option. Could also look into garbage collection.
        % N.B.: DENOISE/ICA_AROMA_component_assessment.pdf might be worth
        % keeping
        
        % done!

		        
    case 'checkrequirements'
        
         
        % sanity checks:
        
        if (~isfield(aap.directory_conventions,'fsldir') ...
            || isempty(aap.directory_conventions.fsldir) ...
                || ~exist(aap.directory_conventions.fsldir,'dir'))
                
                aas_log(aap, true, sprintf('%s: aap.directory_conventions.fsldir is not set properly. Exiting...',mfilename));
        end
        
        % check if AROMA is installed
        [~, AROMA] = aas_cache_get(aap,'aroma');
        AROMA_FNAME = fullfile(AROMA.toolPath,'ICA_AROMA.py');
        if ~exist(AROMA_FNAME,'file')               
            aas_log(aap, true, sprintf('%s: ICA_AROMA must be installed in %s. Exiting...', mfilename, AROMA_FNAME));
        end

        % AROMA requires python3.6
        [~,w] = aas_shell('python --version'); 
        verPy = regexp(w,'(?<=Python )[0-9\.]*','match');
        if ~startsWith(verPy{1},'3.6'), aas_log(aap,true,sprintf('%s: ICA_AROMA requires python 3.6 (not found). Exiting...', mfilename)); end

        % may as well check for template while we're here
        
        WARPREF = fullfile(aap.directory_conventions.fsldir,'data/standard/MNI152_T1_2mm.nii.gz');             

        if ~exist(WARPREF,'file')               
            aas_log(aap, true, sprintf('%s: cannot find FSL template %s. Exiting...', mfilename, WARPREF));
        end
       
    otherwise
        
        aas_log(aap, 1, sprintf('%s: Unknown task %s',mfilename, task));
        
end