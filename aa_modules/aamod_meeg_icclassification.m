function [aap, resp] = aamod_meeg_icclassification(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        aap = aas_report_addimage(aap,subj,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '.jpg']));
        aap = aas_report_add(aap,subj,'<table><tr><th>Accepted</th><th>Rejected</th></tr><tr>');
        for sfx = {'accepted','rejected'}
            aap = aas_report_add(aap,subj,'<td valign="top">');
            for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),['^diagnostic_.*' sfx{1} '.*jpg$']))'
                aap=aas_report_addimage(aap,subj,fn{1});
            end
            aap = aas_report_add(aap,subj,'</td>');
        end
        aap = aas_report_add(aap,subj,'</tr></table>');        
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [junk, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        EEG = pop_loadset(infname{strcmp(spm_file(infname,'ext'),'set')});
        
        switch aas_getsetting(aap,'method')
            case 'ICLabel'
                EEG = iclabel(EEG, 'default');
                prob = aas_getsetting(aap,'criteria.prob');
                goodIcIdx = true(size(EEG.etc.ic_classification.ICLabel.classifications,1),1);
                for p = strsplit(prob,':')
                    % opration
                    if any(p{1}(1) == '*+'), op = p{1}(1); p{1}(1) = []; else, op = '*'; end
                    if op == '*', op = @times; else, op = str2func(op); end
                    
                    % criterion
                    crit = regexp(p{1},'[\>\=\<]'); crit = p{1}(crit);
                    expr = strsplit(p{1},crit);
                    if crit == '=', crit = @eq; else, crit = str2func(crit); end
                    
                    % calculate 'goodness'
                    pIc = crit(EEG.etc.ic_classification.ICLabel.classifications(:,strcmp(EEG.etc.ic_classification.ICLabel.classes,expr{1})),str2double(expr{2}));
                    goodIcIdx = op(goodIcIdx,pIc);
                end
                goodIcIdx  = find(goodIcIdx);
                if ~isempty(isempty(EEG.dipfit)) && ~isempty(aas_getsetting(aap,'criteria.rv'))
                    rvList    = [EEG.dipfit.model.rv];
                    goodRvIdx = find(rvList < aas_getsetting(aap,'criteria.rv'));
                    finalIcIdx = intersect(goodIcIdx, goodRvIdx);
                else
                    finalIcIdx = goodIcIdx;
                end                
            otherwise
                aas_log(aap,true,sprintf('Method %s not yet implemented',aas_getsetting(aap,'method')))
        end
        
        
        pop_viewprops(EEG,0,1:size(EEG.icaweights,1),{},{},0,aas_getsetting(aap,'method'));
        set(gcf,'position',[0,0,1080 1080]);
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-noui',fullfile(aas_getsesspath(aap,subj,sess),sprintf('diagnostic_%s.jpg',mfilename)),'-djpeg','-r150');
        close(gcf);

        
        % Plot rejected components
        sfx = {'rejected' 'accepted'};
        for ic = 1:size(EEG.icaweights,1)
            f = pop_prop_extended(EEG, 0, ic, NaN, {}, {}, 0, aas_getsetting(aap,'method'));
            set(f,'PaperPositionMode','auto');
            print(f,'-noui',fullfile(aas_getsesspath(aap,subj,sess),sprintf('diagnostic_%s_%s_IC%03d.jpg',mfilename,sfx{any(finalIcIdx==ic)+1},ic)),'-djpeg','-r300');
            close(f);
        end
        
        EEG = pop_subcomp(EEG, finalIcIdx, 0, 1);
        EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(finalIcIdx,:);
        
        % save
        outfname = spm_file(infname,'prefix','icclass_');
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',spm_file(outfname{1},'basename'));

        EL.unload;
                
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end