function [aap, resp] = aamod_meeg_fei(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        % for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),'^diagnostic_.*jpg$'))'
        %     aap = aas_report_add(aap,subj,'<table><tr><td>');
        %     aap=aas_report_addimage(aap,subj,fn{1});
        %     aap = aas_report_add(aap,subj,'</td></tr></table>');
        % end
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [~, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        indfnEEG = strcmp(spm_file(infname,'ext'),'set');
        EEG = pop_loadset('filepath',spm_file(infname{indfnEEG},'path'),'filename',spm_file(infname{indfnEEG},'filename'));
        EL.unload;

        data = cell2mat(arrayfun(@(ch) abs(hilbert(EEG.data(ch,:)))', 1:EEG.nbchan, 'UniformOutput',false));

        %% DFA
        cfg = aas_getsetting(aap,'dfa');
        if ~isempty(cfg.windowsize) && ~isempty(cfg.windowoverlap)
            expDFA = calculateDFA(data, cfg.windowsize*EEG.srate, cfg.windowoverlap);
        else
            expDFA = NaN(EEG.nbchan,1);
        end

        %% fEI
        cfg = aas_getsetting(aap,'fei');
        fEI = calculateFEI(data, cfg.windowsize*EEG.srate, cfg.windowoverlap);

        %% Save
        [~, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;

        data = keepfields(ft_read_header(infname{indfnEEG}),{'elec','label'});
        data.dimord = 'chan';

        FT.unload;

        % dfa
        if ~all(isnan(expDFA))
            dfa = data;
            dfa.expdfa = expDFA;
            save(spm_file(infname{indfnEEG},'prefix','dfa_','ext','mat'),'dfa');
            aas_desc_outputs(aap,'meeg_session',[subj sess],'dfa',spm_file(infname{indfnEEG},'prefix','dfa_','ext','mat'));
        end

        % fei
        fei = data;
        fei.fei = fEI;
        save(spm_file(infname{indfnEEG},'prefix','fei_','ext','mat'),'fei');
        aas_desc_outputs(aap,'meeg_session',[subj sess],'fei',spm_file(infname{indfnEEG},'prefix','fei_','ext','mat'));

    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
