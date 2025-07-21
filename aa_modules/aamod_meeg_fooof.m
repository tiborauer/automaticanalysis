function [aap, resp] = aamod_meeg_fooof(aap,task,subj)

resp='';

switch task
    case 'report'

    case 'doit'
        [~, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;

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

            % Save into FT data structure
            load(fnTF{1},'timefreq');
            r2 = load(spm_file(fnTF{1},'prefix','fooof_','suffix','_rsquared','ext',''));

            % - aperiodic
            ap = load(spm_file(fnTF{1},'prefix','fooof_','suffix','_aperiodic','ext',''));
            fooof = keepfields(timefreq,{'label' 'elec'});
            fooof.dimord = 'chan';
            fooof.r2 = cellfun(@(ch) r2.(ch), fooof.label);
            fooof.offset = cellfun(@(ch) ap.(ch)(1), fooof.label);
            fooof.exp = cellfun(@(ch) ap.(ch)(2), fooof.label);            
            save(spm_file(fnTF{1},'prefix','fooof_','suffix','_aperiodic','ext',''),'fooof');

            % - peaks
            peaks = load(spm_file(fnTF{1},'prefix','fooof_','suffix','_peaks','ext',''));
            fooof = keepfields(timefreq,{'label' 'elec'});
            fooof.dimord = 'chan_band';
            fooof.bandspec = aas_getsetting(aap,'bandspec');
            fooof.band = fieldnames(fooof.bandspec);
            fooof.r2 = cellfun(@(ch) r2.(ch), fooof.label);
            fooof.peakfreq = cell(0,0);
            fooof.peakbandwidth = cell(0,0);
            fooof.peakpower = cell(0,0);
            for b = 1:numel(fooof.band)
                for ch = 1:numel(fooof.label)
                    sel1 = find(peaks.(fooof.label{ch})(:,1)>=fooof.bandspec.(fooof.band{b})(1), 1, 'first');
                    sel2 = find(peaks.(fooof.label{ch})(:,1)<fooof.bandspec.(fooof.band{b})(2), 1, 'last');
                    fooof.peakfreq(ch,b) = {peaks.(fooof.label{ch})(sel1:sel2,1)};
                    fooof.peakbandwidth(ch,b) = {peaks.(fooof.label{ch})(sel1:sel2,3)};
                    fooof.peakpower(ch,b) = {peaks.(fooof.label{ch})(sel1:sel2,2)};
                end
            end
            save(spm_file(fnTF{1},'prefix','fooof_','suffix','_peaks','ext',''),'fooof');
        end

        aas_desc_outputs(aap,'subject',subj,'aperiodic',spm_select('FPList',aas_getsubjpath(aap,subj),'^fooof_.*_aperiodic.mat$'));
        aas_desc_outputs(aap,'subject',subj,'peaks',spm_select('FPList',aas_getsubjpath(aap,subj),'^fooof_.*_peaks.mat$'));
        FOOOF.unload;
        FT.unload;
    case 'checkrequirements'
end
end