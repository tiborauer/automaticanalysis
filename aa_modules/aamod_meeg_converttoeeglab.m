function [aap, resp] = aamod_meeg_converttoeeglab(aap,task,subj,sess)

resp='';

switch task
      case 'report'
        % summary
        % - init
        if ~isfield(aap.report, aap.tasklist.currenttask.name)
            aap.report.(aap.tasklist.currenttask.name).pnts = NaN(aas_getN_bydomain(aap,'subject'),numel(aap.acq_details.meeg_sessions));
        end

        % - save pnts
        outFn = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj,sess],'meeg','output'));        
        [~, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        EEG = pop_loadset(outFn{strcmp(spm_file(outFn,'ext'),'set')});
        EL.unload;
        aap.report.(aap.tasklist.currenttask.name).pnts(subj,sess) = EEG.pnts;

        % figures
        for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),'^diagnostic_.*jpg$'))'
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap = aas_report_addimage(aap,subj,fn{1});
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [~, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        % read data
        EEG = [];
        for i = 1:numel(infname)
            try
                EEG = pop_fileio(infname{i});
                res = evalc('eeg_checkset(EEG)');
                if isempty(res), break;
                else, throw(res); end
            catch err
                E(i) = err;
            end
        end
        if isempty(EEG)
            for i = 1:numel(E)
                aas_log(aap,i==numel(E),sprintf('ERROR: reading %s - %s',infname{i}, E(i).message));
            end
        end
        
        % channel layout
        EEG = pop_chanedit(EEG,'lookup',aas_getfiles_bystream(aap,'channellayout'));
        
        % remove channel
        if ~isempty(aas_getsetting(aap,'removechannel'))
            chns = strsplit(aas_getsetting(aap,'removechannel'),':');
            EEG = pop_select(EEG, 'nochannel', cellfun(@(x) find(strcmp({EEG.chanlocs.labels}, x)), chns));
        end
        
        % downsample        
        if ~isempty(aas_getsetting(aap,'downsample',sess)) && ~isnan(aas_getsetting(aap,'downsample',sess))
            sRate = aas_getsetting(aap,'downsample',sess);
            if sRate < EEG.srate, EEG = pop_resample( EEG, sRate); end
        end
        
        % edit
        % - specify operations
        toEditsetting = aas_getsetting(aap,'toEdit');
        toEditsubj = toEditsetting(...
            cellfun(@(x) any(strcmp(x,aas_getsubjname(aap,subj))),{toEditsetting.subject}) | ...
            strcmp({toEditsetting.subject},'*')...
            );        
        toEdit = struct('type',{},'operation',{});
        for s = 1:numel(toEditsubj)
            sessnames = regexp(toEditsubj(s).session,':','split');
            if any(strcmp(sessnames,aas_getsessname(aap,sess))) || strcmp(sessnames{1}, '*')
                toEdit = horzcat(toEdit,toEditsubj(s).event);
            end
        end
        
        % - do it
        if ~isempty(toEdit)
            for e = toEdit
                if ischar(e.type)
                    try ind = ~cellfun(@isempty, regexp({EEG.event.type},e.type)); 
                    catch
                        aas_log(aap,false,sprintf('No trial %s found', e.type));
                        ind = false(size(EEG.event));
                    end
                elseif isnumeric(e.type)
                    ind = e.type;
                end
                op = strsplit(e.operation,':');
                if ~any(ind) && ~startsWith(op{1},'insert'), continue; end
                switch op{1}
                    case 'clear'
                        EEG.event(1:end) = [];
                        EEG.urevent(1:end) = [];
                    case 'remove'
                        urind = [EEG.event(ind).urevent];
                        EEG.event(ind) = [];
                        EEG.urevent(urind) = [];
                    case 'keep'
                        urind = [EEG.event(ind).urevent];
                        EEG.event = EEG.event(ind);
                        EEG.urevent = EEG.urevent(urind);
                    case 'keepbeforeevent'
                        ind = find(ind);
                        indCrit = find(strcmp({EEG.event.type},op{2}));                        
                        indKeep = arrayfun(@(l) find(ind<l, 1, 'last'), indCrit);
                        indOmit = setdiff(1:numel(ind),indKeep);
                        EEG.event(ind(indOmit)) = [];
                        EEG.urevent(ind(indOmit)) = [];
                    case 'rename'
                        for i = find(ind)
                            EEG.event(i).type = op{2};
                            EEG.urevent(i).type = op{2};
                        end
                    case 'unique'
                        ex = [];
                        switch op{2}
                            case 'first'
                                for i = 2:numel(ind)
                                    if ind(i) && ind(i-1), ex(end+1) = i; end
                                end
                            case 'last'
                                for i = 1:numel(ind)-1
                                    if ind(i) && ind(i+1), ex(end+1) = i; end
                                end
                        end
                        EEG.event(ex) = [];
                        EEG.urevent(ex) = [];
                    case 'iterate'
                        ind = cumsum(ind).*ind;
                        for i = find(ind)
                            EEG.event(i).type = sprintf('%s%02d',EEG.event(i).type,ind(i));
                            EEG.urevent(i).type = sprintf('%s%02d',EEG.urevent(i).type,ind(i));
                        end
                    case 'insert'
                        loc = str2num(op{2});
                        newE = EEG.event(loc);
                        for i = 1:numel(newE)
                            newE(i).type = e.type;
                        end
                        events = EEG.event(1:loc(1)-1);
                        for i = 1:numel(loc)-1
                            events = [events newE(i) EEG.event(loc(i):loc(i+1)-1)];
                        end
                        if isempty(i), i = 0; end
                        events = [events newE(i+1) EEG.event(loc(i+1):end)];
                        EEG.event = events;
                        EEG.urevent = rmfield(events,'urevent');
                    case {'insertwithlatency' 'insertwithtime'}
                        latencies = op{2};
                        if ~any(strcmp(latencies,{'beginning' 'end'})), latencies = str2num(latencies); end
                        if ischar(latencies) % special cases
                            switch latencies
                                case 'beginning'
                                    latencies = 1;
                                case 'end'
                                    latencies = EEG.pnts;
                            end
                        else
                            if strcmp(op{1}, 'insertwithtime')
                                % find the samples at latencies exact or just later
                                latencies = arrayfun(@(lat) find(EEG.times - lat>=0, 1, 'first'), latencies);
                            end
                        end
                        newE = struct(...
                            'type',e.type,...
                            'duration',EEG.srate/1000,...
                            'timestamp',[],...
                            'latency',num2cell(latencies),...
                            'urevent',0 ...
                            );
                        if ~isfield(EEG.event,'timestamp'), newE = rmfield(newE,'timestamp'); end
                        EEG.event = [EEG.event newE];
                        [~, ord] = sort([EEG.event.latency]);
                        EEG.event = EEG.event(ord);                        
                        for i = 1:numel(EEG.event), EEG.event(i).urevent = i; end
                        EEG.urevent = rmfield(EEG.event,'urevent');
                    case 'inserteachbetween'
                        ind_start = find(strcmp({EEG.event.type},op{2}),1,'first');
                        ind_end = find(strcmp({EEG.event.type},op{4}),1,'last');
                        latencies = EEG.event(ind_start).latency:str2num(op{3})*EEG.srate:EEG.event(ind_end).latency;
                        latencies(1) = latencies(1)+1; % first sample after start
                        if latencies(end) == EEG.event(ind_end).latency, latencies(end) = latencies(end)-1; end % last sample before end if overlap
                        for lat = latencies
                            EEG.event(end+1) = struct('type',e.type,...
                                                      'duration',EEG.srate/1000,...
                                                      'timestamp',[],...
                                                      'latency',lat,...
                                                      'urevent',[]);
                        end
                        EEG.event(strcmp({EEG.event.type},'empty')) = [];
                        [~, ord] = sort([EEG.event.latency]);
                        EEG.event = EEG.event(ord);                        
                        for i = 1:numel(EEG.event), EEG.event(i).urevent = i; end
                        EEG.urevent = rmfield(EEG.event,'urevent');
                    case 'prefixpattern'
                        pfx = strsplit(op{3}(2:end-1),' ');
                        if sum(ind) ~= numel(pfx), aas_log(aap, true, sprintf('You need %d prefixes in a space-delimited list',sum(ind))); end
                        ind = find(ind);
                        for n = 1:numel(ind)
                            EEG.event(ind(n)).type = regexprep(EEG.event(ind(n)).type,['(' op{2} ')'],[pfx{n} '$1']);
                            EEG.urevent(EEG.event(ind(n)).urevent).type = EEG.event(ind(n)).type;
                        end
                    case 'ignorebefore'
                        EEG = pop_select(EEG,'nopoint',[0 EEG.event(find(ind,1,'first')).latency-1]);
                        beInd = find(strcmp({EEG.event.type},'boundary'),1,'first');
                        samplecorr = EEG.event(beInd).duration;
                        EEG.event(beInd) = [];
                        ureindcorr = EEG.event(1).urevent -1;
                        
                        % adjust events                        
                        for i = 1:numel(EEG.event)
                            EEG.event(i).urevent = EEG.event(i).urevent - ureindcorr;
                        end
                        EEG.urevent(1:ureindcorr) = [];
                        
                        % adjust time
                        for i = 1:numel(EEG.urevent)
                            EEG.urevent(i).latency = EEG.urevent(i).latency - samplecorr;
                        end
                    case 'ignoreafter'
                        if EEG.pnts > EEG.event(find(ind,1,'last')).latency
                            EEG = pop_select(EEG,'nopoint',[EEG.event(find(ind,1,'last')).latency+1 EEG.pnts]);
                            beInd = find(strcmp({EEG.event.type},'boundary'),1,'last');
                            EEG.event(beInd) = [];
                        end
                    otherwise
                        aas_log(aap,false,sprintf('Operation %s not yet implemented',op{1}));
                end
                % update events
                for i = 1:numel(EEG.event)
                    urind = find(strcmp({EEG.urevent.type},EEG.event(i).type) & [EEG.urevent.latency]==EEG.event(i).latency);
                    if isempty(urind) % try based on timing only
                        urind = find([EEG.urevent.latency]==EEG.event(i).latency);urind = find(strcmp({EEG.urevent.type},EEG.event(i).type) & [EEG.urevent.latency]==EEG.event(i).latency);
                    end
                    EEG.event(i).urevent = urind(1);
                end
            end
        end

        % diagnostics
        diagpath = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_raw.jpg']);
        meeg_diagnostics_continuous(EEG,aas_getsetting(aap,'diagnostics'),sprintf('Raw with %d timepoint (~%d s)',EEG.pnts, round(EEG.xmax-EEG.xmin)),diagpath);
                
        fnameroot = sprintf('eeg_%s',aas_getsubjname(aap,subj));
        while ~isempty(spm_select('List',aas_getsesspath(aap,subj,sess),[fnameroot '.*']))
            fnameroot = ['aa' fnameroot];
        end
        
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',fnameroot);
        outfname = spm_select('FPList',aas_getsesspath(aap,subj,sess),[fnameroot '.*']);
        
        EL.unload;
        
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end