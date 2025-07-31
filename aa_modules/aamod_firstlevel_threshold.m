% AA module - first level thresholding
% **********************************************************************
% You should no longer need to change this module - you may just
% modify the .xml or model in your user script
% **********************************************************************
% Tibor Auer MRC CBU Cambridge 2012-2013
%
% CHANGE HISTORY
%
% 07/2018 --  added explicit template parameter to xml. Save SPM
% stats table. Added zero-sig-voxel watermark. Save all maps and 
% renders (even zero sig voxel results). 
% Add sanity check(s). Added optional "description" text which is
% overlayed on map if defined. General cleanup. [MSJ]

function [aap,resp]=aamod_firstlevel_threshold(aap,task,subj)

resp='';

switch task
    
    case 'report'
        
        % collect contrast names and prepare summary
        contrasts = aas_getsetting(aas_setcurrenttask(aap,aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream(1).sourcenumber),'contrasts');
        cons = [contrasts(2:end).con];
        conNames = {cons.name};
        [~,a] = unique(conNames,'first');
        conNames = conNames(sort(a));
        
        if subj == 1 % first
            for C = 1:numel(conNames)
                if  ~isfield(aap.report,sprintf('html_C%02d',C))
                    aap.report.(sprintf('html_C%02d',C)).fname = fullfile(aap.report.condir,[aap.report.fbase sprintf('_C%02d.htm',C)]);
                    aap = aas_report_add(aap,'C00',...
                        sprintf('<a href="%s" target=_top>%s</a><br>',...
                        aap.report.(sprintf('html_C%02d',C)).fname,...
                        ['Contrast: ' conNames{C}]));
                    aap = aas_report_add(aap,sprintf('C%02d',C),['HEAD=Contrast: ' conNames{C}]);
                end
                if ~isempty(aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix)
                    aap = aas_report_add(aap,sprintf('C%02d',C),sprintf('<h2>Branch: %s</h2>',...
                        aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix(2:end)));
                end
            end
        end
        
        fnSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        load(fnSPM,'SPM');

        % sanity check -- make sure SPM.swd has the correct path
        if ~isequal(SPM.swd, spm_file(fnSPM,'path')), SPM.swd = spm_file(fnSPM,'path'); end
        
        for C = 1:numel(SPM.xCon)
            
            conName = strrep_multi(SPM.xCon(C).name,{' ' ':' '>'},{'' '_' '-'});
            conInd = find(strcmp(conNames,SPM.xCon(C).name));
            if isempty(conInd), continue, end
            
            aap = aas_report_add(aap,sprintf('C%02d',conInd),['Subject: ' basename(aas_getsubjpath(aap,subj)) '<br>']);
            
            aap = aas_report_add(aap,subj,sprintf('<h4>%02d. %s</h4>',conInd,conName));
            
            f{1} = fullfile(aas_getsubjpath(aap,subj),...
                sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_overlay_3_001.jpg',conInd,conName));
            
            % older versions didn't create overlay/renders if no voxels
            % survived thresholding, ergo the check here. We now create
            % all images, but this check doesn't hurt, and may be useful
            % if generating a report on an old extant analysis
            
            if exist(f{1},'file')
                
                tstat = dlmread(strrep(f{1},'_overlay_3_001.jpg','.txt'));
                
                f{2} = fullfile(aas_getsubjpath(aap,subj),...
                    sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_render.jpg',conInd,conName));
                
                % add overlay and render images to single subject report...
                
                aap = aas_report_add(aap, subj,'<table><tr>');
                aap = aas_report_add(aap, subj, sprintf('T = %2.2f - %2.2f</tr><tr>', tstat(1), tstat(2)));
                for i = 1:2
                    aap = aas_report_add(aap, subj,'<td>');
                    aap = aas_report_addimage(aap, subj, f{i});
                    aap = aas_report_add(aap, subj,'</td>');
                end
                
                % add SPM stats table
                statsfname = fullfile(aas_getsubjpath(aap,subj),sprintf('table_firstlevel_threshold_C%02d_%s.jpg', conInd, conName));
                if ~exist(statsfname,'file')
                    make_stats_table(SPM, statsfname, C, ...
                        aap.tasklist.currenttask.settings.threshold.p, ...
                        aap.tasklist.currenttask.settings.threshold.correction);
                end
                aap = aas_report_add(aap, subj,'<td>');
                aap = aas_report_addimage(aap, subj, statsfname);
                aap = aas_report_add(aap, subj,'</td>');
                
                aap = aas_report_add(aap,subj,'</tr></table>');
                
                % ...also add images & table to module report
                
                aap = aas_report_add(aap,sprintf('C%02d',conInd),'<table><tr>');
                aap = aas_report_add(aap,sprintf('C%02d',conInd),sprintf('T = %2.2f - %2.2f</tr><tr>', tstat(1), tstat(2)));
                for i = 1:2
                    aap = aas_report_add(aap, sprintf('C%02d',conInd),'<td>');
                    aap = aas_report_addimage(aap,sprintf('C%02d',conInd), f{i});
                    aap = aas_report_add(aap,sprintf('C%02d',conInd),'</td>');
                end
                aap = aas_report_add(aap, sprintf('C%02d',conInd),'<td>');
                aap = aas_report_addimage(aap, sprintf('C%02d',conInd), statsfname);
                aap = aas_report_add(aap,sprintf('C%02d',conInd),'</td>');
                aap = aas_report_add(aap,sprintf('C%02d',conInd),'</tr></table>');
                
            end
            
        end
        
    case 'doit'
        
        % sanity checks
        
        if (strcmp(aap.tasklist.currenttask.settings.overlay.template,'structural'))
            aas_log(aap, false, sprintf('WARNING (%s): You should verify template ''structural'' is in the same space as ''epi''.', mfilename));
        end
        
        % Init
        
        try doTFCE = aap.tasklist.currenttask.settings.threshold.doTFCE; catch, doTFCE = 0; end % TFCE?
        corr = aap.tasklist.currenttask.settings.threshold.correction;		% correction
        u0   = aap.tasklist.currenttask.settings.threshold.p;			% height threshold
        nSl = aap.tasklist.currenttask.settings.overlay.nth_slice;
        tra = aap.tasklist.currenttask.settings.overlay.transparency;
        Outputs.thr = '';
        Outputs.cl = '';
        Outputs.sl = '';
        Outputs.Rend = '';
        
        cwd=pwd;
        localroot = aas_getsubjpath(aap,subj);
        anadir = fullfile(localroot,aap.directory_conventions.stats_singlesubj);
        cd(anadir);
        
        % we now explicitly define which structural template (native, normed, or SPMT1)
        % to use for overlay mapping rather than leaving it up to provenance
        
        switch  aap.tasklist.currenttask.settings.overlay.template
            
            case 'structural'
                
                % we can guarantee normalised_structural and SPMT1, but "structural"
                % could be native or normalised depending on tasklist. That's why
                % the option is called "structural" not "native" (although native
                % is intended).	If the user picks this option, we assume they know
                % what they're doing...
                
                if aas_stream_has_contents(aap,subj,'structural')
                    tmpfile = aas_getfiles_bystream(aap, subj, 'structural');
                else
                    aas_log(aap, true, sprintf('%s: Cannot find structural. Exiting...', mfilename));
                end
                
            case 'SPMT1'
                
                % assume a reasonable default location, but assume the user put
                % the correct location in aap.dir_con.T1template if it's not empty
                
                tmpfile = 'toolbox/OldNorm/T1.nii';
                if ~isempty(aap.directory_conventions.T1template), tmpfile = aap.directory_conventions.T1template; end
                if (tmpfile(1) ~= '/'), tmpfile = fullfile(fileparts(which('spm')),tmpfile); end
                
                if ~exist(tmpfile,'file')
                    aas_log(aap, true, sprintf('%s: SPM T1 template not found. Exiting...', mfilename));
                end
                
            otherwise
                
                aas_log(aap, true, sprintf('%s: Unknown template option. Exiting...', mfilename));
                
        end
        
        % Now get contrasts...
        
        SPM=[];
        fnSPM = aas_getfiles_bystream(aap, subj,'firstlevel_spm');
        load(fnSPM,'SPM');
        
        for c = 1:numel(SPM.xCon)
            no_sig_voxels = false; % need this for later
            conName = strrep_multi(SPM.xCon(c).name,{' ' ':' '>'},{'' '_' '-'});
            STAT = SPM.xCon(c).STAT;
            df = [SPM.xCon(c).eidf SPM.xX.erdf];
            XYZ  = SPM.xVol.XYZ;
            S    = SPM.xVol.S;   % Voxel
            R    = SPM.xVol.R;   % RESEL
            V = spm_vol(fullfile(anadir, SPM.xCon(c).Vspm.fname));
            Z = spm_get_data(SPM.xCon(c).Vspm,XYZ);
            dim = SPM.xCon(c).Vspm.dim;
            VspmSv   = cat(1,SPM.xCon(c).Vspm);
            n = 1; % No conjunction
            
            if doTFCE
                k = aas_getsetting(aap,'threshold.extent');
                job.spmmat = {fnSPM};
                job.mask = {fullfile(fileparts(fnSPM),'mask.nii,1')};
                job.conspec = struct( ...
                    'titlestr','', ...
                    'contrasts',c, ...
                    'n_perm',5000, ...
                    'vFWHM',0 ...
                    );
                job.tbss = 0;
                job.openmp = 1;
                cg_tfce_estimate(job);
                iSPM = SPM;
                iSPM.title = '';
                iSPM.Ic = c;
                iSPM.stattype = 'TFCE';
                iSPM.thresDesc = corr;
                iSPM.u = u0;
                iSPM.k = k;
                [SPM, xSPM] = cg_get_tfce_results(iSPM);
                Z = xSPM.Z;
                XYZ = xSPM.XYZ;
                if isempty(Z)
                    aas_log(aap,false,sprintf('INFO: No voxels survive TFCE(%s)=%1.4f, k=%0.2g',corr, u0, k));
                    no_sig_voxels = true;
                end
            else
                % Height threshold filtering
                switch corr
                    case 'iTT'
                        % TODO
                        [Z, XYZ, th] = spm_uc_iTT(Z,XYZ,u0,1);
                    case 'FWE'
                        u = spm_uc(u0,df,STAT,R,n,S);
                    case 'FDR'
                        u = spm_uc_FDR(u0,df,STAT,n,VspmSv,0);
                    case 'none'
                        u = spm_u(u0^(1/n),df,STAT);
                end
                Q      = find(Z > u);
                Z      = Z(:,Q);
                XYZ    = XYZ(:,Q);
                if isempty(Q)
                    aas_log(aap,false,sprintf('INFO: No voxels survive height threshold u=%0.2g',u));
                    no_sig_voxels = true;
                end
                
                % Extent threshold filtering
                if ischar(aas_getsetting(aap,'threshold.extent')) % probability-based
                    k = strsplit(aas_getsetting(aap,'threshold.extent'),':'); k{2} = str2double(k{2});
                    iSPM = SPM;
                    iSPM.Ic = c;
                    iSPM.thresDesc = corr;
                    iSPM.u = u0;
                    iSPM.k = 0;
                    iSPM.Im = [];
                    [~,xSPM] = spm_getSPM(iSPM);
                    T = spm_list('Table',xSPM);
                    switch k{1}
                        case {'FWE' 'FDR'}
                            k{1} = ['p(' k{1} '-corr)'];
                        case {'none'}
                            k{1} = 'p(unc)';
                    end
                    pInd = strcmp(T.hdr(1,:),'cluster') & strcmp(T.hdr(2,:),k{1});
                    kInd = strcmp(T.hdr(2,:),'equivk');
                    k = min(cell2mat(T.dat(cellfun(@(p) ~isempty(p) && p<k{2}, T.dat(:,pInd)),kInd)));
                    if isempty(k), k = Inf; end
                else
                    k = aas_getsetting(aap,'threshold.extent');
                end
                
                A     = spm_clusters(XYZ);
                Q     = [];
                for i = 1:max(A)
                    j = find(A == i);
                    if length(j) >= k
                        Q = [Q j];
                    end
                end
                Z     = Z(:,Q);
                XYZ   = XYZ(:,Q);
                if isempty(Q)
                    aas_log(aap,false,sprintf('INFO: No voxels survive extent threshold k=%0.2g',k));
                    no_sig_voxels = true;
                end
            end
            
            % Reconstruct
            Yepi  = zeros(dim(1),dim(2),dim(3));
            indx = sub2ind(dim,XYZ(1,:)',XYZ(2,:)',XYZ(3,:)');
            Yepi(indx) = Z;
            V.fname = spm_file(V.fname,'basename',strrep(spm_file(V.fname,'basename'),'spm','thr'));
            V.descrip = sprintf('thr{%s_%1.4f;ext_%d}%s',corr,u0,k,V.descrip(strfind(V.descrip,'}')+1:end));
            spm_write_vol(V,Yepi);
            
            % Cluster
            clusterfname = spm_file(V.fname,'prefix','cl_');
            if ~isempty(aas_getsetting(aap,'cluster'))
                if no_sig_voxels
                    copyfile(V.fname,clusterfname);
                else
                    switch aas_getsetting(aap,'cluster.method')
                        case 'fusionwatershed'
                            [s,FWS] = aas_cache_get(aap,'fws');
                            if ~s
                                aas_log(aap,false,'WARNING: Fusion-Watershed is not installed! --> clustering skipped');
                            else
                                FWS.load;
                                settings = aas_getsetting(aap,'cluster.options.fusionwatershed');
                                obj = fws.generate_ROI(V.fname,...
                                    'threshold_method','z','threshold_value',0.1,...
                                    'filter',settings.extentprethreshold,'radius',settings.searchradius,'merge',settings.mergethreshold,...
                                    'plot',false,'output',true);
                                
                                % - exclude small (<k) ROIs
                                smallROIs = obj.table.ROIid(obj.table.Volume < k);
                                obj.label(reshape(arrayfun(@(l) any(l == smallROIs), obj.label(:)), obj.grid.d)) = 0;
                                obj.table(arrayfun(@(l) any(l == smallROIs), obj.table.ROIid),:) = [];
                                    
                                % - save results
                                save.vol(obj.label,obj.grid,spm_file(clusterfname,'ext',''),'Compressed',false);
                                writetable(obj.table,spm_file(clusterfname,'ext','csv'));
                                FWS.unload;
                            end
                    end
                end
            end          
            
            % Overlay
            % - edges of activation
            slims = ones(4,2);
            sAct = arrayfun(@(x) any(reshape(Yepi(x,:,:),[],1)), 1:size(Yepi,1));
            if numel(find(sAct))<2, slims(1,:) = [1 size(Yepi,1)];
            else, slims(1,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
            sAct = arrayfun(@(y) any(reshape(Yepi(:,y,:),[],1)), 1:size(Yepi,2));
            if numel(find(sAct))<2, slims(2,:) = [1 size(Yepi,2)];
            else, slims(2,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
            sAct = arrayfun(@(z) any(reshape(Yepi(:,:,z),[],1)), 1:size(Yepi,3));        
            if numel(find(sAct))<2, slims(3,:) = [1 size(Yepi,3)];
            else, slims(3,:) = [find(sAct,1,'first') find(sAct,1,'last')]; end
            % - convert to mm
            slims = sort(V.mat*slims,2);
            % - extend if too narrow (min. 50mm)
            slims = slims + (repmat([-25 25],4,1).*repmat(diff(slims,[],2)<50,1,2));         

            % - draw
            axis = {'sagittal','coronal','axial'};
            for a = 1:3
                if ~no_sig_voxels, stat_fname = {V.fname}; else, stat_fname = {}; end
                [fig, v] = map_overlay(tmpfile,stat_fname,axis{a},slims(a,1):nSl:slims(a,2));
                fnsl{a} = fullfile(localroot, sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_overlay_%d.jpg',c,conName,a));
                
                if (~isempty(aap.tasklist.currenttask.settings.description))
                    annotation('textbox',[0 0.5 0.5 0.5],'String',aap.tasklist.currenttask.settings.description,'FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                end
                
                if (no_sig_voxels)
                    annotation('textbox',[0 0.475 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
                end
                
                fnsl{a} = spm_print(fnsl{a},fig,'jpg');
            end
            
            dlmwrite(fullfile(localroot, sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s.txt',c,conName)),[min(v(v~=0)), max(v)]);
            
            % Render
            
            % FYI: render should always work regardless of template type because it
            % maps input into MNI, if necessary
            
            if ~no_sig_voxels
                if numel(Z)  < 2 % render fails with only one active voxel
                    Z = horzcat(Z,Z);
                    XYZ = horzcat(XYZ,XYZ);
                end
                % render fails with single first slice
                for a = 1:3
                    if all(XYZ(a,:)==1)
                        Z = horzcat(Z,Z(end));
                        XYZ = horzcat(XYZ,XYZ(:,end)+circshift([1;0;0],a-1));
                    end
                end
            end
            
            dat.XYZ = XYZ;
            dat.t = Z';
            dat.mat = SPM.xVol.M;
            dat.dim = dim;
            rendfile  = aap.directory_conventions.Render;
            if ~exist(rendfile,'file') && (rendfile(1) ~= '/'), rendfile = fullfile(fileparts(which('spm')),rendfile); end
            fn3d = fullfile(localroot,sprintf('diagnostic_aamod_firstlevel_threshold_C%02d_%s_render.jpg',c,conName));
            global prevrend
            prevrend = struct('rendfile',rendfile, 'brt',0.5, 'col',eye(3));
            out = spm_render(dat,0.5,rendfile); spm_figure('Close','Graphics');
            img = vertcat(horzcat(out{1},out{3},out{5}),horzcat(out{2},out{4},out{6}));
            fig = figure;
            imshow(img,'Border','tight');
            
            if (~isempty(aap.tasklist.currenttask.settings.description))
                annotation('textbox',[0 0.5 0.5 0.5],'String',aap.tasklist.currenttask.settings.description,'FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
            end
            
            if (no_sig_voxels)
                annotation('textbox',[0 0.45 0.5 0.5],'String','No voxels survive threshold','FitBoxToText','on','fontweight','bold','color','y','fontsize',18,'backgroundcolor','k');
            end
            
            print(fig,'-noui',fn3d,'-djpeg','-r300');
            close(fig);
            
            % Outputs
            
            if exist(V.fname,'file'), Outputs.thr = strvcat(Outputs.thr, V.fname); end
            if exist(clusterfname,'file'), Outputs.cl = strvcat(Outputs.cl, clusterfname); end
            for f = 1:numel(fnsl)
                if exist(fnsl{f},'file'), Outputs.sl = strvcat(Outputs.sl, fnsl{f}); end
            end
            if exist(fn3d,'file'), Outputs.Rend = strvcat(Outputs.Rend, fn3d); end
            
        end
        
        cd (cwd);
        
        % Describe outputs
        
        aap=aas_desc_outputs(aap,subj,'firstlevel_thr',Outputs.thr);
        aap=aas_desc_outputs(aap,subj,'firstlevel_clusters',Outputs.cl);
        aap=aas_desc_outputs(aap,subj,'firstlevel_thrslice',Outputs.sl);
        aap=aas_desc_outputs(aap,subj,'firstlevel_thr3D',Outputs.Rend);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 0, sprintf('%s: Unknown task %s',mfilename, task));
end

end