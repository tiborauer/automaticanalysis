function [aap resp]=aamod_diffusion_dtinlfit(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        %% Fetch inputs
        % Get nii filenames from stream
        diffinput=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data');
        bvals=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvals');
        bvecs=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvecs');
        betmask=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'BETmask');
        
        % Find which line of betmask contains the brain mask
        for betind=1:size(betmask,1)
            if strfind(betmask(betind,:),'bet_nodif_brain_mask')
                break
            end;
        end;        
        betmask=betmask(betind,:);
        
        %% Apply dtinlfit
        data_in = spm_read_vols(spm_vol(diffinput));
        data_mask = spm_read_vols(spm_vol(betmask));
        bval = dlmread(bvals);
        bvec = dlmread(bvecs);
    
        if strcmp(aap.options.wheretoprocess, 'localsingle') && ~isempty(aap.directory_conventions.poolprofile) % local execution - use parfor for slices
            profiles = parallel.clusterProfiles;
            if ~any(strcmp(profiles,aap.directory_conventions.poolprofile))
                ppfname = which(spm_file(aap.directory_conventions.poolprofile,'ext','.settings'));
                if isempty(ppfname)
                    aas_log(aap,true,sprintf('ERROR: settings for pool profile %s not found!',aap.directory_conventions.poolprofile));
                else
                    P=parcluster(parallel.importProfile(ppfname));
                end
            else
                aas_log(aap,false,sprintf('INFO: pool profile %s found',aap.directory_conventions.poolprofile));
                P=parcluster(aap.directory_conventions.poolprofile);
            end
            switch class(P)
                case 'parallel.cluster.Torque'
                    aas_log(aap,false,'INFO: Torque engine is detected');
                    P.ResourceTemplate = '-l nodes=^N^,mem=100MB,walltime=00:05:00';
                case 'parallel.cluster.Local'
                    aap.options.aaparallel.numberofworkers = min([12 feature('numCores')]);
            end
            aas_makedir(aap, fullfile(aas_getsesspath(aap,subjind,diffsessind),'cluster'));
            P.JobStorageLocation = fullfile(aas_getsesspath(aap,subjind,diffsessind),'cluster');
            P.NumWorkers = min([...
                size(data_mask,3) ... % for each slice
                aap.options.aaparallel.numberofworkers
                ]);

            % submit
            for z = 1:size(data_in,3)
                job(z) = createJob(P);
                if isprop(job(z),'AutoAttachFiles'), job(z).AutoAttachFiles = false; end
                createTask(job(z), @dti_slice, 7,{squeeze(data_in(:,:,z,:)), data_mask(:,:,z), bval, bvec});
                submit(job(z));
            end
            % monitor
            while ~all(strcmp({job.State},'finished') | strcmp({job.State},'failed')), pause(1); end
            if any(strcmp({job.State},'failed')), aas_log(aap,true,['Slice(s) ' num2str(find(strcmp({job.State},'failed'))) ' failed']); end
            %retrieve
            for z = 1:size(data_in,3)
                out = fetchOutputs(job(z));
                [S0(:,:,z), L1(:,:,z), L2(:,:,z), L3(:,:,z), V1(:,:,z,:), V2(:,:,z,:), V3(:,:,z,:)] = out{:};
            end
        else % within-module cluster computing is not available
            for z = 1:size(data_in,3)
                [S0(:,:,z), L1(:,:,z), L2(:,:,z), L3(:,:,z), V1(:,:,z,:), V2(:,:,z,:), V3(:,:,z,:)] = dti_slice(squeeze(data_in(:,:,z,:)), data_mask(:,:,z), bval, bvec);
            end
        end
        
        dti.S0 = S0;
        dti.L1 = L1;
        dti.L2 = L2;
        dti.L3 = L3;
        dti.V1 = V1;
        dti.V2 = V2;
        dti.V3 = V3;
        
        dti.AD = L1;
        dti.RD = (L2+L3)/2;
        dti.MD = (L1+L2+L3)/3;
        dti.FA = sqrt(3/2)*...
            sqrt(var(cat(4, L1, L2, L3),[],4)*2./...
            (L1.^2+L2.^2+L3.^2));
        
        % Now describe outputs
        V = spm_vol(betmask); V.dt = spm_type('float32');
        sesspath = aas_getpath_bydomain(aap,'diffusion_session',[subjind,diffsessind]);

        outstreams=aas_getstreams(aap,'output');        
        for outind=1:length(outstreams)
            metric = strrep(outstreams{outind},'dti_','');
            if ~isfield(dti,metric)
                aas_log(aap,false,sprintf('Metric %s for stream %s not exist!',metric,outstreams{outind}));
                continue; 
            end
            Y = dti.(metric);  
            nifti_write(fullfile(sesspath,[outstreams{outind} '.nii']),Y,outstreams{outind},V);
            aap=aas_desc_outputs(aap,'diffusion_session',[subjind,diffsessind],outstreams{outind},[outstreams{outind} '.nii']);
        end
end
end

%% DTI wrapper

function [S0, L1, L2, L3, V1, V2, V3] = dti_slice(data, mask, b, bv)

S0 = mask*0;
L1 = mask*0;
L2 = mask*0;
L3 = mask*0;
V1 = repmat(mask*0,[1 1 3]);
V2 = repmat(mask*0,[1 1 3]);
V3 = repmat(mask*0,[1 1 3]);

for x = 1:size(data,1)
    for y = 1:size(data,2)
        
        vdata = squeeze(data(x,y,:))';
        
        if ~mask(x,y), continue; end
        
        [S0(x,y), DT] = dti_nlfit(vdata, b, bv);
        
        [V,D] = eig(DT);

        L1(x,y) = D(3,3);
        L2(x,y) = D(2,2);
        L3(x,y) = D(1,1);
        
        V1(x,y,:) = V(:,3);
        V2(x,y,:) = V(:,2);
        V3(x,y,:) = V(:,1);
    end
end
end