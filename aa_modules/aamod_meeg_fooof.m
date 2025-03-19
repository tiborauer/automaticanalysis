function [aap, resp] = aamod_meeg_fooof(aap,task,subj)

resp='';

switch task
    case 'report'

    case 'doit'
        % Config fooof
        [~, FOOOF] = aas_cache_get(aap,'fooof');
        FOOOF.load;
        condasetup = ['source ' aap.directory_conventions.condasetup '; conda activate ' FOOOF.condaEnvironment];
        py_script = which(FOOOF.script);
        py_args = sprintf('--freq_range %1.3f %1.3f --aperiodic_mode %s', ...
            aas_getsetting(aap,'frequencyrange'), ...
            aas_getsetting(aap,'aperiodicmode'));
        
        % Run fooof
        for fnTF = cellstr(aas_getfiles_bystream(aap,'subject',subj,'timefreq'))'
            [~, w] = aas_shell([condasetup ';python ' py_script ' ' py_args ' '  fnTF{1} ' ' spm_file(fnTF{1},'prefix','fooof_','ext','')]);
            aas_log(aap,false,w);
        end
        
        aas_desc_outputs(aap,'subject',subj,'aperiodic',spm_select('FPList',aas_getsubjpath(aap,subj),'^fooof_.*_aperiodic.mat$'));
        aas_desc_outputs(aap,'subject',subj,'peaks',spm_select('FPList',aas_getsubjpath(aap,subj),'^fooof_.*_peaks.mat$'));
        FOOOF.unload;
    case 'checkrequirements'
end
end