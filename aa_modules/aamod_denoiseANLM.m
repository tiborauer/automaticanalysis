% AA module
% Denoise a series of structurals...
% NOTE: This is based upon and requires the MRI Denoising Package by Pierrick Coupe
% (This is not provided with AA at this moment)
% Please cite his work, as described within the MRIDenoisingAONLM function!

function [aap,resp]=aamod_denoiseANLM(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 align';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        Sfn = aas_getfiles_bystream(aap,subj,'structural');
        
        try
            DCMfn = aas_getfiles_bystream(aap,subj,'structural_dicom_header');
            
            % dcmhdr{n}.SeriesDescription
            dcmhdr = [];
            load(DCMfn);
        catch
        end
        
        %% Denoise the images
        
        outstruct = '';
        outresid = '';
        
        for d = aap.tasklist.currenttask.settings.structurals
            aas_log(aap,false,sprintf('Denoise structural image %s!', Sfn(d,:)))
            try
                aas_log(aap,false,sprintf('\t structural type %s!', dcmhdr{d}.SeriesDescription))
            catch
                aas_log(aap,false,sprintf('\t structural type UNKNOWN!'))
            end
            
            V = spm_vol(Sfn(d,:));
            Y = spm_read_vols(V);
            
            % Denoised Y
            dY = MRIDenoisingAONLM(Y, ...
                aap.tasklist.currenttask.settings.patchsize, ...
                aap.tasklist.currenttask.settings.searcharea, ...
                aap.tasklist.currenttask.settings.beta, ...
                aap.tasklist.currenttask.settings.rician, ...
                aap.tasklist.currenttask.settings.verbose);
            
            % Residuals...
            rY = dY - Y;
            
            % Save filtered image & residual image
            dV = V;
            [pth, fn, ext] = fileparts(dV.fname);
            dV.fname = fullfile(pth, ['d' fn ext]);
            outstruct = strvcat(outstruct, dV.fname);
            spm_write_vol(dV, dY);
            
            rV = V;
            [pth, fn, ext] = fileparts(rV.fname);
            rV.fname = fullfile(pth, ['res_d' fn ext]);
            outresid = strvcat(outresid, rV.fname);
            spm_write_vol(rV, rY);
            
            try close(2); catch; end
            % Then plot a diagnostic image...
            % Save graphical output to common diagnostics directory
            if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
                mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
            end

            try
                %% Draw noisy and denoised structural...
                spm_check_registration(strvcat(V.fname, dV.fname))
                
                spm_orthviews('reposition', [0 0 0])
                
                figure(1);
                set(gcf,'PaperPositionMode','auto')
                print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' aap.acq_details.subjects(subj).subjname '_' num2str(d) '.jpeg']));
            catch
            end
            if aap.tasklist.currenttask.settings.diagnostic
                
                movieFilename = fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' aap.acq_details.subjects(subj).subjname '_' num2str(d) '.avi']);
                % Create movie file by defining aviObject
                try delete(movieFilename); catch; end
                aviObject = avifile(movieFilename,'compression','none');
                
                try close(2); catch; end
                figure(2)
                set(2, 'Position', [0 0 1000 600])
                windowSize = get(2,'Position');
                
                % Get resliced structural
                EPIlims = [min(Y(:)) max(Y(:))];
                
                for z = 1:size(Y,3)
                    h = subplot(1,2,1);
                    imagesc(rot90(squeeze(Y(:,:,z))))
                    caxis(EPIlims)
                    axis equal off
                    title('Before denoiseANLM')
                    zoomSubplot(h, 1.2)
                    
                    h = subplot(1,2,2);
                    imagesc(rot90(squeeze(dY(:,:,z))))
                    caxis(EPIlims)
                    axis equal off
                    title('After denoiseANLM')
                    zoomSubplot(h, 1.2)
                    
                    % Capture frame and store in aviObject
                    pause(0.01)
                    aviObject = addframe(aviObject,getframe(2,windowSize));
                end
                
                junk = close(aviObject);
            end
        end
        
        %% DESCRIBE OUTPUTS!
        
        try
            dcmhdr = {dcmhdr{aap.tasklist.currenttask.settings.structurals}};
            save(DCMfn, 'dcmhdr')
            
            aap=aas_desc_outputs(aap,subj,'structural_dicom_header', DCMfn);
        catch
        end
        
        % Structural image after denoising
        aap=aas_desc_outputs(aap,subj,'structural', outstruct);
        
        % Residual image after denoising
        aap=aas_desc_outputs(aap,subj,'denoiseResidual', outresid);
        
end