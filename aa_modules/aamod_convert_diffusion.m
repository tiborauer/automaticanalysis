% Get diffusion dicom and convert nii.gz
% Converting Data to Analyze format and extracting the Gradient Directions
% extract the gradient direction and b-values as a text file

function [aap resp]=aamod_convert_diffusion(aap,task,subjind,diffsessind,phaseencodeind)
resp='';

switch task
    case 'report'
    case 'doit'
        
        if exist('phaseencodeind','var')
            indices = [subjind diffsessind phaseencodeind];
            domain = 'diffusion_session_phaseencode_direction';
        else
            indices = [subjind diffsessind];
            domain = 'diffusion_session';
        end;
        domainpath=aas_getpath_bydomain(aap,domain,indices);
        
        % Get DICOM filenames from stream and bvals and bvecs from the
        % header
        [aap niifiles DICOMHEADERS subdirs]=aas_convertseries_fromstream(aap,domain,indices,'dicom_diffusion');
        
        % Now move dummy scans to dummy_scans directory
        
        % From header of this module
        if isfield(aap.tasklist.currenttask.settings,'numdummies') &&...
                ~isempty(aap.tasklist.currenttask.settings.numdummies)
            ndummies=aap.tasklist.currenttask.settings.numdummies;
        else % backward compatibility
            ndummies = aap.acq_details.numdummies;
        end
        
        dummylist=[];
        if ndummies
            dummypath=fullfile(domainpath,'dummy_scans');
            aap=aas_makedir(aap,dummypath);
            for dummyind=1:ndummies
                cmd=['mv ' niifiles{dummyind} ' ' dummypath];
                [pth nme ext]=fileparts(niifiles{dummyind});
                dummylist=strvcat(dummylist,fullfile('dummy_scans',[nme ext]));
                [s w]=aas_shell(cmd);
                if (s)
                    aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n',niifiles{dummyind},dummypath));
                end
            end
        end
        niifiles = {niifiles{ndummies+1:end}};
        DICOMHEADERS = {DICOMHEADERS{ndummies+1:end}};
        % 4D conversion [TA]
        for fileind=1:numel(niifiles)
            V(fileind)=spm_vol(niifiles{fileind});
        end
        if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D
            niifiles = niifiles{1};
            ind = find(niifiles=='-');
            niifiles = [niifiles(1:ind(2)-1) '.nii'];
			spm_file_merge(char({V.fname}),niifiles,0);
        end
        
        % bvals and bvecs based on Guy William's algorithm as implemented by Matthew Brett
        % Init
        y_flipper = diag([1 -1 1]);
        
        % Get voxel to dicom rotation matrix
        orient           = reshape(DICOMHEADERS{1}.ImageOrientationPatient,[3 2]);
        orient(:,3)      = null(orient');
        if det(orient)<0, orient(:,3) = -orient(:,3); end;
        vox_to_dicom = orient;
        vox_to_dicom = vox_to_dicom * y_flipper;
        
        n_hdrs = numel(DICOMHEADERS);
        bvals = zeros(n_hdrs, 1);
        bvecs = zeros(n_hdrs, 3);

        for h = 1:n_hdrs
            % Read B_matrix info
            bm = aas_get_numaris4_numval(DICOMHEADERS{h}.CSAImageHeaderInfo,'B_matrix')';
            % If no B_matrix, this is 0 B value
            if isempty(bm)
                continue
            end
            B_matrix = [bm(1:3); bm(2) bm(4) bm(5); bm(3) bm(5) bm(6)];
            % find max eigenvalue, eigenvector from B_matrix
            [vecs, vals] = eig(B_matrix);
            vals = max(vals);
            [bvals(h), i] = max(vals);
            dbvec = vecs(:,i);
            % For convenience, turn vectors to point towards positive X
            if dbvec(1) < 0
                dbvec = dbvec * -1;
            end
            bvecs(h,:) = [vox_to_dicom\dbvec]';
        end
        
%         bvecs=zeros(length(dicomheader),3);
%         for headerind=1:length(dicomheader)
%                 bvals(headerind)=aas_get_numaris4_numval(dicomheader{headerind}.CSAImageHeaderInfo,'B_value');
%                 if bvals(headerind)~=0
%                     bvecs(headerind,:)=aas_get_numaris4_numval(dicomheader{headerind}.CSAImageHeaderInfo,'DiffusionGradientDirection');
%                 end;
%         end;
        
        % Output final data
        sesspth=aas_getpath_bydomain(aap,domain,indices);
        
        % Write bvals
        bvals_fn=fullfile(sesspth,'bvals');
        fid=fopen(bvals_fn,'w');
        fprintf(fid,'%d ',bvals);
        fprintf(fid,'\n');
        fclose(fid);
        
        % Write bvecs
        bvecs_fn=fullfile(sesspth,'bvecs');
        fid=fopen(bvecs_fn,'w');
        for ln=1:3
            fprintf(fid,'%.14f ',bvecs(:,ln));
            fprintf(fid,'\n');
        end;
        fclose(fid);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,domain,indices,'dummyscans',dummylist);
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_data',niifiles);
        dcmhdrfn = fullfile(domainpath,'dicom_headers.mat');
        save(dcmhdrfn,'DICOMHEADERS');
        aap=aas_desc_outputs(aap,domain,indices,'diffusion_dicom_header',dcmhdrfn);
        aap=aas_desc_outputs(aap,domain,indices,'bvals',bvals_fn);
        aap=aas_desc_outputs(aap,domain,indices,'bvecs',bvecs_fn);    
       
end
end

