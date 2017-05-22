% Automated Defacing Tools (Freesurfer software)
% It defaces structural (T1) and produces a mask which can be applied to any coregistered image
%
% N.B.: It requires FreeSurfer templates talairach_mixed_with_skull.gca and
% face.gca downloaded from https://surfer.nmr.mgh.harvard.edu/fswiki/mri_deface,
% exracted and stored in aap.directory_conventions.templatedir as
% freesurfer_deface_talairach_mixed_with_skull.gca and
% freesurfer_deface_face.gca, respectively.

function [aap resp]=aamod_freesurfer_deface(aap,task,subjind)
resp='';

switch task
    case 'report'
    case 'doit'
        % Check templates
        global aa
        if ~isa(aa,'aaClass'), aa = aaClass; end
        
        tmp_face = fullfile(aap.directory_conventions.templatedir,'freesurfer_deface_face.gca');
        tmp_skull = fullfile(aap.directory_conventions.templatedir,'freesurfer_deface_talairach_mixed_with_skull.gca');
        if ~exist(tmp_face,'file') || ~exist(tmp_skull,'file')
            aas_log(aap,true,sprintf('Templates required: %s, %s\n',tmp_face, tmp_skull));
        end
        
        % Get input
        sdir = pwd;
        Simg = aas_getfiles_bystream(aap,subjind,'structural'); 
        [p, fn, ext] = fileparts(Simg);
        out = fullfile(p,['defaced_' fn ext]);
        
        % Apply mri_deface
        cd(p);
        cmd = sprintf('mri_deface %s %s %s %s',Simg,tmp_skull,tmp_face,out);
        [s, w]=aas_runFScommand(aap,cmd);
        if (s)
            aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
        end;
        cd(sdir);
        
        % Create mask (for other images)
        inf = spm_vol(out);
        Y = spm_read_vols(inf);
        Y = Y>0;
        out_mask = fullfile(p,['defaced_mask_' fn ext]);
        nifti_write(out_mask,Y,'Defaced mask',inf)
        
        % Now describe outputs
        aap=aas_desc_outputs(aap,subjind,'defaced_structural',out);
        aap=aas_desc_outputs(aap,subjind,'defaced_mask',out_mask);
    case 'checkrequirements'
        tmp_face = fullfile(aap.directory_conventions.templatedir,'freesurfer_deface_face.gca');
        tmp_skull = fullfile(aap.directory_conventions.templatedir,'freesurfer_deface_talairach_mixed_with_skull.gca');
        if ~exist(tmp_face,'file') || ~exist(tmp_skull,'file')
            msg = ['ERROR: module requires FreeSurfer templates freesurfer_deface_talairach_mixed_with_skull.gca and freesurfer_deface_face.gca\n'...
                '    1. Download <a href = "https://surfer.nmr.mgh.harvard.edu/pub/dist/mri_deface/talairach_mixed_with_skull.gca.gz">talairach_mixed_with_skull.gca.gz</a> and <a href = "https://surfer.nmr.mgh.harvard.edu/pub/dist/mri_deface/face.gca.gz">face.gca.gz</a>\n'...
                '    2. Extract them\n'...
                '    3. Rename them to freesurfer_deface_talairach_mixed_with_skull.gca and freesurfer_deface_face.gca, respectively\n'...
                sprintf('    4. Move them to %s',aap.directory_conventions.templatedir)];
            aas_log(aap,false,msg);
        end
end
end



