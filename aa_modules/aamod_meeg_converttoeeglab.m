function [aap, resp] = aamod_meeg_converttoeeglab(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),'^diagnostic_.*jpg$'))'
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fn{1});
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [junk, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        % read data
        EEG = [];
        for i = 1:numel(infname)
            try
                EEG = pop_fileio(infname{i});
                [res,o1,o2] = evalc('eeg_checkset(EEG)');
                if isempty(res), break;
                else, throw(res); end
            catch err
                E(i) = err;
            end
        end
        if isempty(EEG)
            for i = 1:numel(E)
                aas_log(aap,false,sprintf('ERROR: reading %s - %s',infname{i}, E(i).message));
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
        if ~isempty(aas_getsetting(aap,'downsample'))
            sRate = aas_getsetting(aap,'downsample');
            if sRate ~= EEG.srate, EEG = pop_resample( EEG, aas_getsetting(aap,'downsample')); end
        end
        
        
        
        % behaviour
        % - grab subject specific variables from aap structure
        behaviourSetting = aas_getsetting(aap,'behaviour'); % as cell array
        behaviourSubj = behaviourSetting(strcmp({behaviourSetting.subject},aas_getsubjname(aap,subj)));
        if isempty(behaviourSubj)
            behaviourSubj = behaviourSetting(strcmp({behaviourSetting.subject},'*')); 
        end    
        
        behaviourSess = behaviourSubj(strcmp({behaviourSubj.session},aas_getsessname(aap,sess)));
        if isempty(behaviourSess)
            behaviourSess = behaviourSubj(strcmp({behaviourSubj.session},'*')); 
        end  
        
        % if, else
        if numel(behaviourSess)>1
            aas_log(aap,false,sprintf('ERROR: task structure for %s',aas_getsessdesc(aap,subj,sess))); % false = warning, true = error
        elseif isempty(behaviourSess)
            % error
        else
            structtask = behaviourSess.taskstructure;

            % - do it
            for j=1:size(structtask,2)
                structtask(j).latency = structtask(j).picstart * 1000;
            end

            % Find all events with event marker of stimulus
            allEventTypes = {EEG.event.type}'; %'
            eventIdx = find(strcmp(allEventTypes, behaviourSess.eventtype));

            % I and J to loop through event marker indexes and table indexes
            I=zeros(1,numel(eventIdx))'; %' preallocate
            for k=1:numel(eventIdx)
                I(k) = k;
            end
            J = eventIdx;

            % Integrate information into each stimulus marker
            for k=1:numel(I)
                i = I(k); % absolute number of epoch event
                j = J(k); % index number of epoch event

                EEG.event(j).picnum = structtask(i).picnum;
                EEG.event(j).soundnum = structtask(i).soundnum;
                EEG.event(j).trialType = structtask(i).trialType;
                EEG.event(j).correct = structtask(i).correct;
                EEG.event(j).RT = structtask(i).RTs;
            end
        end


        
        % edit
        % - specify operations
        toEditsetting = aas_getsetting(aap,'toEdit'); % cell array, 
        toEditsubj = toEditsetting(...
            cellfun(@(x) any(strcmp(x,aas_getsubjname(aap,subj))),{toEditsetting.subject}) | ...
            strcmp({toEditsetting.subject},'*')...
            );        
        toEdit = struct('type',{},'operation',{});
        for s = 1:numel(toEditsubj)
            sessnames = regexp(toEditsubj(s).session,':','split');
            if any(strcmp(sessnames,aas_getsessname(aap,sess))) || sessnames{1} == '*'
                toEdit = horzcat(toEdit,toEditsubj(s).event);
            end
        end
        
        % - do it
        if ~isempty(toEdit)
            for e = toEdit
                if ischar(e.type)
                    ind = ~cellfun(@isempty, regexp({EEG.event.type},e.type));
                elseif isnumeric(e.type)
                    ind = e.type;
                end
                op = strsplit(e.operation,':');
                if ~any(ind) && ~strcmp(op{1},'insert'), continue; end
                switch op{1}
                    case 'remove'
                        EEG.event(ind) = [];
                        EEG.urevent(ind) = [];
                    case 'keep'
                        EEG.event = EEG.event(ind);
                        EEG.urevent = EEG.urevent(ind);
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
                    otherwise
                        aas_log(aap,false,sprintf('Operation %s not yet implemented',op{1}));
                end
            end
            for i = 1:numel(EEG.event)
                EEG.event(i).urevent = i;
            end
        end

        % diagnostics
        diagpath = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_raw.jpg']);
        meeg_diagnostics_continuous(EEG,aas_getsetting(aap,'diagnostics'),'Raw',diagpath);
                
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