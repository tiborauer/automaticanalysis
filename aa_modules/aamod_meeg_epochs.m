function [aap, resp] = aamod_meeg_epochs(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        [~, EL] = aas_cache_get(aap,'eeglab');
        
        MAXNTRIAL = 1000; % do not expect more than 1000 trials
        
        outfname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg','output'));
        outfname = outfname(strcmp(spm_file(outfname,'ext'),'set'));
        segments = reshape(str2double(unique(regexp(spm_file(outfname,'basename'),'(?<=seg-)[0-9]+','match','once'))),1,[]);        
        conds = regexp(spm_file(outfname(endsWith(spm_file(outfname,'basename'),'seg-1')),'basename'),'(?<=_)[A-Z-0-9]+','match','once');
        condnum = cellfun(@(x) str2double(regexp(x,'(?<=-)[0-9]+','match')), conds);
        
        % init summary
        % - first session
        if ~isfield(aap.report, mfilename)
            aap.report.(mfilename).alltrials = cell(aas_getN_bydomain(aap,'subject'),numel(aap.acq_details.meeg_sessions));
            aap.report.(mfilename).condcount = cell(aas_getN_bydomain(aap,'subject'),numel(aap.acq_details.meeg_sessions));
        end
        
        %% Individual
        % collect trials
        switch EL.status
            case 'defined', EL.load;
            case 'unloaded', EL.reload;
        end
        for c = 1:numel(conds)
            subjtrials{c} = logical(sparse(numel(segments),MAXNTRIAL));
            for s = segments   
                segoutfname = outfname(endsWith(spm_file(outfname,'basename'),sprintf('seg-%d',s)));
                if numel(segoutfname) ~= numel(conds), aas_log(aap,true,'All segments MUST have all conditions'); end

                eeg = pop_loadset(segoutfname{c});
                
                trials = struct2table(rmfield(eeg.event,intersect(fieldnames(eeg.event),{'duration' 'timestamp' 'latency' 'epoch'})));
                trials.type = cellfun(@get_eventvalue, cellstr(trials.type));
                trials(trials.type ~= condnum(c),:) = [];

                utrials = cellfun(@get_eventvalue, {eeg.urevent.type})';
                indtrials = cumsum(utrials == condnum(c));

                subjtrials{c}(s,indtrials([trials.urevent])) = true;                
            end
            subjtrials{c}(:,max(indtrials)+1:end) = [];
        end
        EL.unload;

        % plot
        fn = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_' aap.acq_details.meeg_sessions(sess).name '_trials.jpg']);
        if ~exist(fn,'file')
            f = figure;
            f.Position = [0 0 200*numel(conds) 200*numel(segments)];
            tiledlayout(1,numel(conds),'TileSpacing','none')
            for c = 1:numel(conds)
                ax = nexttile;
                hold on
                arrayfun(@(s) plot(find(subjtrials{c}(s,:)),-s,'b*','MarkerSize',7), 1:numel(segments))
                hold off
                ylim([-max(segments+0.5),-0.5]); yticks(-max(segments):-1);
                if c == 1, yticklabels(arrayfun(@(x) sprintf('segment %d',x),max(segments):-1:1,'UniformOutput',false)); 
                else, yticklabels({});
                end
                xlim([0 size(subjtrials{c},2)]);
                title(ax, conds{c});
            end            
            print(f,'-djpeg','-r300',fn);
            close(f);
        end
        aap=aas_report_addimage(aap,subj,fn);
        aap.report.(mfilename).alltrials{subj,sess} = subjtrials;
        
        % table
        aap = aas_report_add(aap,subj,'<table id="data"><tr>');
        aap = aas_report_add(aap,subj,'<th>Segment</th>');
        for c = 1:numel(conds)
            aap = aas_report_add(aap,subj,sprintf('<th>%s</th>',conds{c}));
        end
        aap = aas_report_add(aap,subj,'</tr>');
        for s = segments
            aap = aas_report_add(aap,subj,'<tr>');
            aap = aas_report_add(aap,subj,sprintf('<td>segment %d</td>',s));
            for c = 1:numel(conds)
                condcount(s,c) = nnz(subjtrials{c}(s,:));
                aap = aas_report_add(aap,subj,sprintf('<td>%d</td>',condcount(s,c)));
            end
            aap = aas_report_add(aap,subj,'</tr>');
        end
        aap = aas_report_add(aap,subj,'</table>');
        condcount(condcount==0) = NaN; % zero epoch should be the ones omitted
        aap.report.(mfilename).condcount{subj,sess} = condcount;
               
    case 'summary'        
            % missing data
            isTrialMissing = cellfun(@isempty, aap.report.(mfilename).alltrials(:,sess));
            isCondcountMissing = cellfun(@(cc) isempty(cc) || (isscalar(cc) && isnan(cc)), aap.report.(mfilename).condcount(:,sess));
            lastSubjCondcount = find(~isCondcountMissing,1,'last');

            % update conditions and sessions based on the last subject with no missing data
            outfname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[lastSubjCondcount sess],'meeg','output'));
            outfname = outfname(strcmp(spm_file(outfname,'ext'),'set'));
            segments = reshape(str2double(unique(regexp(spm_file(outfname,'basename'),'(?<=seg-)[0-9]+','match','once'))),1,[]);
            conds = regexp(spm_file(outfname(endsWith(spm_file(outfname,'basename'),'seg-1')),'basename'),'(?<=_)[A-Z-0-9]+','match','once');            
            
            if ~isfield(aap.report.(mfilename),'summarysessions')                                
                [~, aap.report.(mfilename).summarysessions] = aas_getN_bydomain(aap,aap.tasklist.currenttask.domain,lastSubjCondcount);
                stagerepname = aap.tasklist.currenttask.name;
                if ~isempty(aap.tasklist.currenttask.extraparameters)
                    stagerepname = [stagerepname aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix];
                end
                aap = aas_report_add(aap,'er',['<h2>Stage: ' stagerepname '</h2>']);
                aap = aas_report_add(aap,'er','<table><tr>');
            end
            aap.report.(mfilename).summarysessions = setdiff(aap.report.(mfilename).summarysessions,sess);
            
            aap = aas_report_add(aap,'er','<td valign="top">');
            aap = aas_report_add(aap,'er',['<h3>Session: ' aap.acq_details.meeg_sessions(sess).name '</h3>']);
            
            % Boxplot for each condition
            jitter = 0.1; % jitter around position
            jitter = (...
                1+(rand([sum(~isCondcountMissing),size(aap.report.(mfilename).condcount{lastSubjCondcount,sess},1)])-0.5) .* ...
                repmat(jitter*2./[1:size(aap.report.(mfilename).condcount{lastSubjCondcount,sess},1)],sum(~isCondcountMissing),1)...
                ) .* ...
                repmat([1:size(aap.report.(mfilename).condcount{lastSubjCondcount,sess},1)],sum(~isCondcountMissing),1);

            condcountFn = fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_' aap.acq_details.meeg_sessions(sess).name '_conditioncount.jpg']);
            condcountFig = figure; condcountFig.Position = [0 0 200*numel(conds) 600*numel(segments)];            
            tiledlayout(condcountFig,1,numel(conds),'TileSpacing','tight');

            trialcountFn = fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_' aap.acq_details.meeg_sessions(sess).name '_trialcount.jpg']);
            trialcountFig = figure; trialcountFig.Position = [0 0 1080 100*numel(conds)]; 
            tiledlayout(trialcountFig,numel(conds),1,'TileSpacing','tight');  

            for c = 1:size(aap.report.(mfilename).condcount{lastSubjCondcount,sess},2)
                figure(condcountFig); ax = nexttile; hold on;
                boxplot(cell2mat(cellfun(@(cc) cc(:,c), aap.report.(mfilename).condcount(~isCondcountMissing,sess), 'UniformOutput', false)')',...
                    'label',arrayfun(@(x) sprintf('Segment %d',x), 1:size(aap.report.(mfilename).condcount,3), 'UniformOutput', false));
                for seg = 1:size(aap.report.(mfilename).condcount{lastSubjCondcount,sess},1)
                    scatter(jitter(:,seg),cellfun(@(cc) cc(seg,c), aap.report.(mfilename).condcount(~isCondcountMissing,sess)),'k','filled','MarkerFaceAlpha',0.4);
                end
                boxValPlot{c} = getappdata(getappdata(gca,'boxplothandle'),'boxvalplot');
                title(ax,conds{c});

                figure(trialcountFig); ax = nexttile;
                tmp = cellfun(@(trl) full(trl{c}), aap.report.(mfilename).alltrials(~isTrialMissing,sess), 'UniformOutput', false);
                maxNTrials = max(cellfun(@numel, tmp));
                for su = 1:numel(tmp)
                    tmp{su}(end+1:maxNTrials) = false;
                end
                trl = sum(squeeze(shiftdim(cat(3,tmp{:}),2)));
                minTRL = min(trl(:));
                im = image(trl-minTRL+1);
                im.Parent.FontSize = 12;
                yticks(1:size(trl,1)); yticklabels(arrayfun(@(x) sprintf('segment %d',x),1:size(trl,1),'UniformOutput',false));
                colormap(jet(numel(tmp)-minTRL+1));
                cb = colorbar('eastoutside'); set(cb,'XTick',[minTRL numel(tmp)]-minTRL+1.5); set(cb,'XTickLabel',get(cb,'XTick')+minTRL-1.5);
                title(ax,conds{c});
            end

            print(condcountFig,'-djpeg','-r300',condcountFn);
            close(condcountFig);
            aap=aas_report_addimage(aap,'er',condcountFn);

            print(trialcountFig,'-djpeg','-r150',trialcountFn);
            close(trialcountFig);
            aap=aas_report_addimage(aap,'er',trialcountFn);
            
            % Stat table
            aap = aas_report_add(aap,'er','<table id="data"><tr>');
            aap = aas_report_add(aap,'er','<th>Segment</th>');
            for c = 1:numel(conds) % for each condition
                aap = aas_report_add(aap,'er',sprintf('<th>%s [median (IQR)]</th>',conds{c}));
                aap = aas_report_add(aap,'er',sprintf('<th>Outliers</th>'));
            end
            aap = aas_report_add(aap,'er','</tr>');
            for seg = 1:size(aap.report.(mfilename).condcount{subj,sess},1)
                aap = aas_report_add(aap,'er','<tr>');
                aap = aas_report_add(aap,'er',sprintf('<td>segment %d</td>',seg));
                for c = 1:numel(conds) % for each condition
                    aap = aas_report_add(aap,'er',sprintf('<td>%3.3f (%3.3f)</td>',boxValPlot{c}.q2(seg),boxValPlot{c}.q3(seg)-boxValPlot{c}.q1(seg)));
                    subjstr = ' None';
                    if boxValPlot{c}.numFiniteLoOutliers(seg)
                        subjstr = strjoin({aap.acq_details.subjects(sort(cell2mat(boxValPlot{c}.outlierrows(seg)))).subjname},' ');
                    end
                    aap = aas_report_add(aap,'er',sprintf('<td>%s</td>',subjstr));
                end
                aap = aas_report_add(aap,'er','</tr>');
            end
            aap = aas_report_add(aap,'er','</table>');
            aap = aas_report_add(aap,'er','</td>');
                        
            if isempty(aap.report.(mfilename).summarysessions), aap = aas_report_add(aap,'er','</tr></table>'); end

    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [~, SPMtool] = aas_cache_get(aap,'spm');
        [~, EL] = aas_cache_get(aap,'eeglab');
        for tbx = {SPMtool EL}
            switch tbx{1}.status
                case 'defined', tbx{1}.load;
                case 'unloaded', tbx{1}.reload;
            end                
        end
        
        EEGLABFILE = infname{strcmp(spm_file(infname,'ext'),'set')};
        EEG = pop_loadset('filepath',spm_file(EEGLABFILE,'path'),'filename',spm_file(EEGLABFILE,'filename'));
       
        % save original timevector
        if isfield(EEG.etc,'clean_sample_mask')
            urtime = nan(1,length(EEG.etc.clean_sample_mask));
            urtime(EEG.etc.clean_sample_mask) = 1;
        else
            urtime = ones(1,EEG.pnts);
        end
        
        % Conditions
        conditions = aas_getsetting(aap,'condition');
        subjmatches=strcmp(aap.acq_details.subjects(subj).subjname,{conditions.subject});
        sessmatches=strcmp(aap.acq_details.meeg_sessions(sess).name,{conditions.session});
        % If no exact spec found, try session wildcard, then subject
        % wildcard, then wildcard for both
        if (~any(sessmatches & subjmatches))
            sesswild=strcmp('*',{conditions.session});
            if (any(sesswild & subjmatches))
                sessmatches=sesswild;
            else
                subjwild=strcmp('*',{conditions.subject});
                if (any(sessmatches & subjwild))
                    subjmatches=subjwild;
                else
                    subjmatches=subjwild;
                    sessmatches=sesswild;
                end
            end
        end
        
        % Should now have just one model spec
        conditionnum = sessmatches & subjmatches;
        eventdef = conditions(conditionnum).event;
        
        % determine if there is any segment specification
        segmentdef = [];
        indSegmentDefiniftion = cell_index({eventdef.conditionlabel},'segment');
        if indSegmentDefiniftion(1)
            segmentdef = eventdef(indSegmentDefiniftion);
            eventdef(indSegmentDefiniftion) = [];
        end
        
        % data rejection, correct events if rejection happens together with an event
        eventRejection = aas_getsetting(aap,'rejectionevent');
        
        % generate segemnt definition for the whole dataset if none exists
        if isempty(segmentdef)
            segmentdef.eventvalue = 'begin:end';
            segmentdef.trlshift = [2 -2];
        end
        
        datafn = {};
        for seg = 1:numel(segmentdef)
            
            % get trials
            if isnumeric(segmentdef(seg).eventvalue) % segment range defined as sample number
                segRangeEvenInd = [...
                    find([EEG.event.latency] <= sum(~isnan(urtime(1:segmentdef(seg).eventvalue(1)))),1,'last') ...
                    find([EEG.event.latency] <= sum(~isnan(urtime(1:segmentdef(seg).eventvalue(2)))),1,'last') ...
                    ];
            else % segment range defined as event type
                segevs = strsplit(segmentdef(seg).eventvalue,':');
                % reinsert events if they have been rejected
                for csegevs = segevs(~(strcmp(segevs,'begin') | strcmp(segevs,'end')))
                    if ~any(strcmp({EEG.event.type},csegevs{1})), EEG = reinsert_urevent(EEG,eventRejection,csegevs{1}); end
                end
                if strcmp(segevs{1},'begin'), segRangeEvenInd(1) = 1;
                else, segRangeEvenInd(1) = find(strcmp({EEG.event.type},segevs{1})); end
                if strcmp(segevs{2},'end'), segRangeEvenInd(2) = numel(EEG.event);
                else, segRangeEvenInd(2) = find(strcmp({EEG.event.type},segevs{2})); end
            end
            
            % segment latency
            segShift = segmentdef(seg).trlshift/1000*EEG.srate;
            
            % segment data
            segEEG = pop_select(EEG,'point',[EEG.event(segRangeEvenInd(1)).latency+segShift(1) EEG.event(segRangeEvenInd(2)).latency+segShift(2)]);
            if isempty(EEG.event(segRangeEvenInd(1)).urevent), segRangeEvenInd(1) = find(~strcmp({EEG.event.type},eventRejection),1,'first'); end
            if isempty(EEG.event(segRangeEvenInd(2)).urevent), segRangeEvenInd(2) = find(~strcmp({EEG.event.type},eventRejection),1,'last'); end
            segEEG.urevent = segEEG.urevent(EEG.event(segRangeEvenInd(1)).urevent:EEG.event(segRangeEvenInd(2)).urevent);
            for e = 1:numel(segEEG.event)
                segEEG.event(e).urevent = segEEG.event(e).urevent - EEG.event(segRangeEvenInd(1)).urevent + 1;
            end
            
            % epoch
            if isempty(eventdef)
                epochEEG = segEEG;
                datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),sprintf('epochs_seg-%d.set',seg));
                datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),sprintf('epochs_seg-%d.fdt',seg));
                pop_saveset(epochEEG,datafn{end});
            else
                for ev = eventdef
                    try
                        epochEEG = pop_epoch(segEEG,{ev.eventvalue},(ev.eventwindow + ev.trlshift)/1000);
                    catch E
                        if ~isempty(regexp(E.message,'dataset .* is empty','once'))
                            fclose(fopen(fullfile(aas_getsesspath(aap,subj,sess),'empty'),'w'));
                            continue;
                        else
                            rethrow(E);
                        end
                    end

                    % baseline correction
                    if ~isempty(ev.baselinewindow)
                        if isnumeric(ev.baselinewindow)
                            epochEEG = pop_rmbase(epochEEG,ev.baselinewindow + ev.trlshift);
                        else
                            switch ev.baselinewindow
                                case 'all'
                                    epochEEG = pop_rmbase(epochEEG,[]);
                                otherwise
                                    % NYI
                            end
                        end
                    end

                    % Inster event-of-interest outside eventwindow as event with 0 latency
                    if aas_getsetting(aap,'ensurefieldtripcompatibility')
                        missingEp = arrayfun(@(ep) isempty(ep.eventtype) || ~any(strcmp(cellstr(ep.eventtype), ev.eventvalue)), epochEEG.epoch);
                        if any(missingEp)
                            durEp = round((epochEEG.xmax-epochEEG.xmin)*epochEEG.srate+1);
                            % find valid epochs
                            urSel = find(strcmp({epochEEG.urevent.type}, ev.eventvalue));
                            urLat = [epochEEG.urevent(urSel).latency];
                            urWin = arrayfun(@(ul) ul+ev.eventwindow+ev.trlshift ,urLat ,'UniformOutput' ,false);
                            urLat = [epochEEG.urevent(urSel(cellfun(@(uw) all(epochEEG.etc.clean_sample_mask(uw(1):uw(2))), urWin))).latency];
                            if numel(urLat) > epochEEG.trials
                                aas_log(aap,false,'Excessive events');
                                urInd = 0;
                                for ep = 1:epochEEG.trials
                                    if isempty(epochEEG.epoch(ep).eventurevent)
                                        urInd(end+1) = urInd(end)+1;
                                        continue
                                    end
                                    if iscell(epochEEG.epoch(ep).eventurevent)
                                        epWin = ev.eventwindow +...
                                            epochEEG.urevent(epochEEG.epoch(ep).eventurevent{1}).latency -...
                                            epochEEG.epoch(ep).eventlatency{1};
                                        urInd(end+1) = find(arrayfun(@(lat) (epWin(1) <= lat) && (epWin(2) >= lat), urLat));
                                    else
                                        epWin = ev.eventwindow +...
                                            epochEEG.urevent(epochEEG.epoch(ep).eventurevent).latency -...
                                            epochEEG.epoch(ep).eventlatency;
                                        urInd(end+1) = find(arrayfun(@(lat) (epWin(1) <= lat) && (epWin(2) >= lat), urLat));
                                    end
                                end
                                urInd(1) = [];
                                urLat = urLat(urInd);
                            elseif numel(urLat) < epochEEG.trials
                                aas_log(aap,true,'Missing events');
                            end

                            urLat = urLat(missingEp);
                            [~, urInd] = intersect([epochEEG.urevent.latency], urLat);
                            evLat = (0:durEp:(epochEEG.trials-1)*durEp)-epochEEG.times(1)+1;
                            evLat = evLat(missingEp);

                            % update events
                            evInsert = struct(...
                                'type', ev.eventvalue,...
                                'duration', 1,...
                                'timestamp', [],...
                                'latency', num2cell(evLat),...
                                'urevent', num2cell(reshape(urInd,size(urLat))),...
                                'epoch', num2cell(1:epochEEG.trials)...
                                );
                            evInsert = rmfield(evInsert,setdiff(fieldnames(evInsert),fieldnames(epochEEG.event)));
                            [~, ~, sortInd] = intersect({'epoch' 'latency'},fieldnames(evInsert));
                            epochEEG.event = table2struct(sortrows(struct2table([evInsert epochEEG.event]),sortInd))';

                            % update epochs
                            for ep = 1:epochEEG.trials
                                indE = find([epochEEG.event.epoch] == ep);
                                if isscalar(indE)
                                    epochEEG.epoch(ep) = struct(...
                                        'event',indE,...
                                        'eventtype',epochEEG.event(indE).type,...
                                        'eventduration',1,...
                                        'eventtimestamp',[],...
                                        'eventlatency',epochEEG.event(indE).latency+epochEEG.xmin*epochEEG.srate-1-(ep-1)*durEp,...
                                        'eventurevent',epochEEG.event(indE).urevent...
                                        );
                                else
                                    epochEEG.epoch(ep).event = indE;
                                    epochEEG.epoch(ep).eventtype= {epochEEG.event(indE).type};
                                    epochEEG.epoch(ep).eventduration = repmat({1},1,numel(indE));
                                    epochEEG.epoch(ep).eventtimestamp = repmat({[]},1,numel(indE));
                                    epochEEG.epoch(ep).eventlatency = num2cell([epochEEG.event(indE).latency]+epochEEG.xmin*epochEEG.srate-1-(ep-1)*durEp);
                                    epochEEG.epoch(ep).eventurevent = {epochEEG.event(indE).urevent};
                                end
                            end                            
                        end
                    end
                    
                    condfn = sprintf('%s-%d',ev.conditionlabel,get_eventvalue(ev.eventvalue));
                    datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),sprintf('epochs_%s_seg-%d.fdt',condfn,seg));
                    datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),sprintf('epochs_%s_seg-%d.set',condfn,seg));                    
                    pop_saveset(epochEEG,datafn{end});                    
                end
                
                % event matching
                eventmatch = aas_getsetting(aap,'eventmatch');
                eminfname = datafn(endsWith(datafn,sprintf('seg-%d.set',seg)));
                for em = eventmatch
                
                    % eventvalues
                    eventvalues = cellfun(@(ev) eventdef(strcmp({eventdef.conditionlabel},ev)).eventvalue, em{1}, 'UniformOutput', false);
                    
                    % data
                    alleeg = cellfun(@(t) pop_loadset(eminfname{contains(eminfname, t)}), em{1});
                
                    if numel(alleeg) ~= numel(eventvalues), aas_log(aap,true,'All EEGs are expected to have a single condition'); end
                    
                    for indeeg = 1:numel(alleeg)
                        % events
                        selEv = [];
                        evur = [alleeg(indeeg).event(strcmp({alleeg(indeeg).event.type}, eventvalues{indeeg})).urevent];
                        urevur = find(strcmp({alleeg(indeeg).urevent.type},eventvalues{indeeg}));
                        for otherind = setdiff(1:numel(alleeg),indeeg)
                            othereegeventsel = strcmp({alleeg(otherind).event.type},eventvalues{otherind});
                            othereegur = [alleeg(otherind).event(othereegeventsel).urevent];
                            for ur = othereegur
                                if strcmp(eventvalues{indeeg}, eventvalues{otherind})
                                    selEv = [selEv find(evur == ur)];
                                elseif otherind < indeeg
                                    try selEv = [selEv find(evur == urevur(find(urevur > ur, 1, 'first')))]; catch, end
                                elseif otherind > indeeg
                                    try selEv = [selEv find(evur == urevur(find(urevur < ur, 1, 'last')))]; catch, end
                                else
                                    aas_log(aap,true,'Conditions not expected');
                                end
                            end
                        end
                        % - select event matched with all other EEGs
                        selEOI = find(strcmp({alleeg(indeeg).event.type}, eventvalues{indeeg}));
                        selEv = tabulate(selEv);
                        selEv = selEv(selEv(:,2) == (numel(alleeg)-1),1);
                
                        % epochs
                        selEp = arrayfun(@(e)...
                            find(arrayfun(@(ep)...
                                any(getIfCellMat(ep.eventurevent)==alleeg(indeeg).event(selEOI(e)).urevent),...
                                alleeg(indeeg).epoch)),...
                            selEv); % assume one event-of-interest per epoch
                
                        % remove from EEG
                        % - adjust event latencies
                        for indSelEOI = selEv'
                            alleeg(indeeg).event(selEOI(indSelEOI)).latency = alleeg(indeeg).event(selEOI(indSelEOI)).latency-...
                                (indSelEOI-sum(selEv<indSelEOI)-1)*...
                                round((alleeg(indeeg).xmax-alleeg(indeeg).xmin)*alleeg(indeeg).srate+1);
                        end
                        alleeg(indeeg).event = alleeg(indeeg).event(selEOI(selEv));
                        alleeg(indeeg).epoch = alleeg(indeeg).epoch(selEp);
                        alleeg(indeeg).data = alleeg(indeeg).data(:,:,selEp);

                        % re-index events with epochs
                        for indE = 1:numel(alleeg(indeeg).event)                            
                            indEp = find(arrayfun(@(ep) any(getIfCellMat(ep.eventurevent) == alleeg(indeeg).event(indE).urevent), alleeg(indeeg).epoch));
                            alleeg(indeeg).event(indE).epoch = indEp;
                            indEvinEp = strcmp(cellstr(alleeg(indeeg).epoch(indEp).eventtype), eventvalues{indeeg});
                            currEp = struct(...
                                'event', indE,...
                                'eventtype', alleeg(indeeg).event(indE).type,...
                                'eventduration',1,...
                                'eventtimestamp',[],...
                                'eventlatency',getIfCellMat(alleeg(indeeg).epoch(indEp).eventlatency,indEvinEp),...
                                'eventurevent',alleeg(indeeg).event(indE).urevent...
                                );
                            alleeg(indeeg).epoch(indEp) = rmfield(currEp,setdiff(fieldnames(currEp),fieldnames(alleeg(indeeg).epoch)));
                        end
                    end
                
                    % save EEGs;
                    for eeg = alleeg
                        datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),['em-' strjoin(em{1},'-') '_' eeg.datfile]);
                        datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),['em-' strjoin(em{1},'-') '_' eeg.filename]);
                        pop_saveset(eeg,'filename',['em-' strjoin(em{1},'-') '_' eeg.filename],'filepath',eeg.filepath);
                    end
                end                
            end               
        end
        
        EL.unload;
        %% Describe outputs
        if any(startsWith(spm_file(datafn,'basename'),'em')) % save only eventmatched
            datafn = datafn(startsWith(spm_file(datafn,'basename'),'em'));
        end
        aap = aas_desc_outputs(aap,subj,sess,'meeg',datafn);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end

function val = get_eventvalue(eventvalue)
if iscell(eventvalue)
    val = cellfun(@(x) get_eventvalue(x), eventvalue);
    return
end

% based on ft_trialfun_general
if isnumeric(eventvalue)
    val = eventvalue;
elseif ischar(eventvalue) && numel(eventvalue)>1 && (eventvalue(1)=='S'|| eventvalue(1)=='R')
    val = str2double(eventvalue(2:end));
else
    warning('Eventvalue %s cannot be parsed',eventvalue);
    val = nan;
end
end

function EEG = reinsert_urevent(EEG,eventRejection,urEvent)
indUre = find(strcmp({EEG.urevent.type},urEvent));
ure = EEG.urevent(indUre);
urtime = nan(1,length(EEG.etc.clean_sample_mask));
urtime(EEG.etc.clean_sample_mask) = 1;

nRej = 0;
ure.latency = round(ure.latency);
if ure.latency == 1
    e = EEG.event(1);
    e.type = ure.type;
    e.latency = 1;
    e.duration = 0;
    e.urevent = indUre;
    EEG.event = [e EEG.event];
else
    if isnan(urtime(ure.latency)) && ~isnan(urtime(ure.latency-1)), nRej = nRej + 1; end
    
    urtime(ure.latency) = 2;
    
    for t = 1:ure.latency
        if isnan(urtime(t))
            if t == 1 || ~isnan(urtime(t-1)), nRej = nRej + 1; end
        end
    end
    
    indRej = find(strcmp({EEG.event.type},eventRejection));
    corrOffset = sum(isnan(urtime(1:ure.latency)));
    e = EEG.event(indRej(nRej));
    e.type = ure.type;
    e.latency = ure.latency - corrOffset;
    e.duration = 0;
    e.urevent = indUre;
    EEG.event = [EEG.event(1:indRej(nRej)) e EEG.event(indRej(nRej)+1:end)];
end
end

function val = getIfCellMat(val, ind)
    if iscell(val), val = cell2mat(val); end
    if nargin > 1, val = val(ind); end
end