function [aap, resp] = aamod_meeg_epochs(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        [~, EL] = aas_cache_get(aap,'eeglab');
        [~, FT] = aas_cache_get(aap,'fieldtrip');

        MAXNTRIAL=1000; % do not expect more than 1000 trials
        
        outfname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg','output'));
        outfname = outfname(strcmp(spm_file(outfname,'ext'),'set'));
        conds = strsplit(spm_file(outfname{1},'basename'),'_'); conds = conds(2:end-1); % first is 'epochs', last is 'seg-x'
        condnum = cellfun(@(x) str2double(regexp(x,'(?<=-)[0-9]+','match')), conds);
        
        % init summary
        % - first session
        if ~isfield(aap.report, mfilename)
            aap.report.(mfilename).trl{1} = cell(aas_getN_bydomain(aap,'subject'),numel(aap.acq_details.meeg_sessions)); % combined trialnumber
            aap.report.(mfilename).trl{2} = cell(aas_getN_bydomain(aap,'subject'),numel(aap.acq_details.meeg_sessions)); % specific trialnumber
            aap.report.(mfilename).condcount = cell(aas_getN_bydomain(aap,'subject'),numel(aap.acq_details.meeg_sessions));
        end
        % - new subject/session
        aap.report.(mfilename).trl{1}{subj,sess} = zeros(numel(outfname),0,numel(conds));
        aap.report.(mfilename).trl{2}{subj,sess} = zeros(numel(outfname),0,numel(conds));
        aap.report.(mfilename).condcount{subj,sess} = nan(numel(outfname),numel(conds));
        
        %% Individual
        % plot
        sfxPlot = {'specific' 'combined'};
        switch EL.status
            case 'defined', EL.load;
            case 'unloaded', EL.reload;
        end 
        for d = 1:numel(outfname)
            EEG(d) = pop_loadset('filepath',spm_file(outfname{d},'path'),'filename',spm_file(outfname{d},'filename'));
        end
        EL.unload;
        switch FT.status
            case 'defined', FT.load;
            case 'unloaded', FT.reload;
        end 
        for d = 1:numel(outfname)
            % remove unwanted events because they may confound the conversion
            EOI = arrayfun(@(x) any(condnum==x), cellfun(@(x) get_eventvalue(x), {EEG(d).event.type}));
            EEG(d).event(~EOI) = [];
            seg(d) = eeglab2fieldtripER(EEG(d));
        end
        FT.unload;        

        for t = 1:2 % specific and combined trialnumber
            trl = nan(numel(outfname),MAXNTRIAL,numel(conds));
            for d = 1:numel(outfname)
                for c = 1:numel(conds)
                    trl(d,table2array(seg(d).ureventinfo(seg(d).trialinfo == condnum(c),t)),c) = -d;
                end
            end
            isTrial = squeeze(sum(isnan(trl),1))<size(trl,1);
            if find(size(isTrial) == MAXNTRIAL) == 2, isTrial = isTrial'; end % make sure it is column
            ntrl = find(sum(isTrial,2),1,'last');
            trl = trl(:,1:ntrl,:);
            fn = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_' aap.acq_details.meeg_sessions(sess).name '_trialcount_' sfxPlot{t} '.jpg']);
            if ~exist(fn,'file')
                f = figure;
                for c = 1:size(trl,3)
                    ax = subplot(1,size(trl,3),c);
                    plot(1:size(trl,2),trl(:,:,c),'*','MarkerSize',7);
                    ylim([-5,0]); yticks(-size(trl,1):-1); yticklabels(arrayfun(@(x) sprintf('segment %d',x),size(trl,1):-1:1,'UniformOutput',false));
                    xlim([0 size(trl,2)+1]); xticks(5:5:size(trl,2)); set(ax,'XGrid','on');
                    title(ax,conds{c});
                end
                
                print(f,'-djpeg','-r300',fn);
                close(f);
            end
            aap=aas_report_addimage(aap,subj,fn);
            aap.report.(mfilename).trl{t}{subj,sess} = ~isnan(trl);
        end
        
        % table
        aap = aas_report_add(aap,subj,'<table id="data"><tr>');
        aap = aas_report_add(aap,subj,'<th>Segment</th>');
        for c = 1:numel(conds)
            aap = aas_report_add(aap,subj,sprintf('<th>%s</th>',conds{c}));
        end
        aap = aas_report_add(aap,subj,'</tr>');
        for d = 1:numel(outfname)
            aap = aas_report_add(aap,subj,'<tr>');
            aap = aas_report_add(aap,subj,sprintf('<td>%s</td>',spm_file(outfname{d},'basename')));
            for c = 1:numel(conds)
                condcount(d,c) = sum(aap.report.aamod_meeg_epochs.trl{1}{subj,sess}(d,:,c));
                aap = aas_report_add(aap,subj,sprintf('<td>%d</td>',condcount(d,c)));
            end
            aap = aas_report_add(aap,subj,'</tr>');
        end
        aap = aas_report_add(aap,subj,'</table>');
        condcount(condcount==0) = NaN; % zero epoch should be the ones omitted
        aap.report.(mfilename).condcount{subj,sess}(:,:) = condcount;
        
        %% Summary in case of more subjects
        if (subj > 1) && (subj == numel(aap.acq_details.subjects)) % last subject
            % missing data
            % - trl
            for t = 1:2
                % -- data size based on full data
                dnd = cellfun(@ndims, aap.report.(mfilename).trl{t}(:,sess));
                sel = dnd == max(dnd);
                ds = max(cell2mat(cellfun(@size, aap.report.(mfilename).trl{t}(sel,sess), 'UniformOutput',false)));
                for dsubj = 1:size(aap.report.(mfilename).trl{t},1)
                    trl = aap.report.(mfilename).trl{t}{dsubj,sess};
                    aap.report.(mfilename).trl{t}{dsubj,sess} = zeros(ds);
                    aap.report.(mfilename).trl{t}{dsubj,sess}(1:size(trl,1),1:size(trl,2),1:size(trl,3)) = trl;
                end
            end

            % - condcount
            isMissing = cellfun(@(cc) isempty(cc) || (isscalar(cc) && isnan(cc)), aap.report.(mfilename).condcount(:,sess));
            ds = size(aap.report.(mfilename).condcount{find(~isMissing,1,'first'),sess});
            aap.report.(mfilename).condcount(isMissing,sess) = {nan(ds)};

            [~, iSess] = aas_getN_bydomain(aap,aap.tasklist.currenttask.domain,subj);
            firstSess = iSess(1);
            lastSess = iSess(end);
            
            if sess == firstSess
                stagerepname = aap.tasklist.currenttask.name;
                if ~isempty(aap.tasklist.currenttask.extraparameters)
                    stagerepname = [stagerepname aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix];
                end
                aap = aas_report_add(aap,'er',['<h2>Stage: ' stagerepname '</h2>']);
                aap = aas_report_add(aap,'er','<table><tr>');
            end
            
            aap = aas_report_add(aap,'er','<td valign="top">');
            aap = aas_report_add(aap,'er',['<h3>Session: ' aap.acq_details.meeg_sessions(sess).name '</h3>']);
            
            % Boxplot for each condition
            jitter = 0.1; % jitter around position
            jitter = (...
                1+(rand([size(aap.report.(mfilename).condcount,1),size(aap.report.(mfilename).condcount{subj,sess},1)])-0.5) .* ...
                repmat(jitter*2./[1:size(aap.report.(mfilename).condcount{subj,sess},1)],size(aap.report.(mfilename).condcount,1),1)...
                ) .* ...
                repmat([1:size(aap.report.(mfilename).condcount{subj,sess},1)],size(aap.report.(mfilename).condcount,1),1);
            condcountFn = fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_' aap.acq_details.meeg_sessions(sess).name '_conditioncount.jpg']);
            condcountFig = figure;
            trialcountFn = fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_' aap.acq_details.meeg_sessions(sess).name '_trialcount.jpg']);
            trialcountFig = figure;
            for c = 1:size(aap.report.(mfilename).condcount{subj,sess},2)
                figure(condcountFig); ax = subplot(1,size(aap.report.(mfilename).condcount{subj,sess},2),c); hold on;
                boxplot(cell2mat(cellfun(@(cc) cc(:,c), aap.report.(mfilename).condcount(:,sess), 'UniformOutput', false)')',...
                    'label',arrayfun(@(x) sprintf('Segment %d',x), 1:size(aap.report.(mfilename).condcount,3), 'UniformOutput', false));
                for s = 1:size(aap.report.(mfilename).condcount{subj,sess},1)
                    scatter(jitter(:,s),cellfun(@(cc) cc(s,c), aap.report.(mfilename).condcount(:,sess)),'k','filled','MarkerFaceAlpha',0.4);
                end
                boxValPlot{c} = getappdata(getappdata(gca,'boxplothandle'),'boxvalplot');
                title(ax,['# condition: ' conds{c}]);
                
                figure(trialcountFig);
                for t = 1:2 % combined and specific trialnumber
                    ax = subplot(2,size(aap.report.(mfilename).condcount{subj,sess},2),(t-1)*size(aap.report.(mfilename).condcount{subj,sess},2)+c);
                    tmp = cellfun(@(trl) trl(:,:,c), aap.report.(mfilename).trl{t}(:,sess), 'UniformOutput', false);
                    trl = squeeze(shiftdim(cat(3,tmp{:}),2));
                    rangeTRL = [min(trl(:)) max(trl(:))];
                    im = image(trl-rangeTRL(1)+1);
                    im.Parent.FontSize = 12;
                    daspect([1 1 1])
                    yticks(1:size(trl,1)); yticklabels(arrayfun(@(x) sprintf('segment %d',x),1:size(trl,1),'UniformOutput',false));
                    colormap(jet(diff(rangeTRL)+1));
                    cb = colorbar('southoutside'); set(cb,'XTick',[rangeTRL(1) median(rangeTRL) rangeTRL(2)]-rangeTRL(1)+1.5); set(cb,'XTickLabel',get(cb,'XTick')+rangeTRL(1)-1.5); cb.Label.String = '# overlapping subjects';
                    title(ax,['# trial: ' conds{c}]);
                end
            end
            set(condcountFig,'Renderer','zbuffer');
            print(condcountFig,'-djpeg','-r300',condcountFn);
            close(condcountFig);
            aap=aas_report_addimage(aap,'er',condcountFn);
            set(trialcountFig,'Renderer','zbuffer');
            print(trialcountFig,'-djpeg','-r150',trialcountFn);
            close(trialcountFig);
            aap=aas_report_addimage(aap,'er',trialcountFn);
            
            % Stat table
            aap = aas_report_add(aap,'er','<table id="data"><tr>');
            aap = aas_report_add(aap,'er','<th>Segment</th>');
            for c = 1:size(aap.report.(mfilename).condcount{subj,sess},2) % for each condition
                if ~isempty(strfind(conds{c},'seg-')), continue; end
                aap = aas_report_add(aap,'er',sprintf('<th>%s [median (IQR)]</th>',conds{c}));
                aap = aas_report_add(aap,'er',sprintf('<th>Outliers</th>'));
            end
            aap = aas_report_add(aap,'er','</tr>');
            for d = 1:size(aap.report.(mfilename).condcount{subj,sess},1)
                aap = aas_report_add(aap,'er','<tr>');
                aap = aas_report_add(aap,'er',sprintf('<td>%s</td>',spm_file(outfname{d},'basename')));
                for c = 1:size(aap.report.(mfilename).condcount{subj,sess},2) % for each condition
                    if ~isempty(strfind(conds{c},'seg-')), continue; end
                    aap = aas_report_add(aap,'er',sprintf('<td>%3.3f (%3.3f)</td>',boxValPlot{c}(d,:).q2,boxValPlot{c}(d,:).q3-boxValPlot{c}(d,:).q1));
                    subjstr = ' None';
                    if boxValPlot{c}(d,:).numFiniteLoOutliers
                        subjstr = strjoin({aap.acq_details.subjects(sort(cell2mat(boxValPlot{c}(d,:).outlierrows))).subjname},' ');
                    end
                    aap = aas_report_add(aap,'er',sprintf('<td>%s</td>',subjstr));
                end
                aap = aas_report_add(aap,'er','</tr>');
            end
            aap = aas_report_add(aap,'er','</table>');
            aap = aas_report_add(aap,'er','</td>');
            
            if sess == lastSess, aap = aas_report_add(aap,'er','</tr></table>'); end
        end

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
        
        % process events
        condfn = '';
        for e = 1:numel(eventdef)
            % detect trialnumber
            condfn = sprintf('%s_%s-%d',condfn,eventdef(e).conditionlabel,get_eventvalue(eventdef(e).eventvalue));
            if e == 1
                % adjust event window
                eventdef(e).eventwindow = (eventdef(e).eventwindow + eventdef(e).trlshift)/1000; % in seconds
                % set baseline correction wondow
                if ~isempty(eventdef(e).baselinewindow) && isnumeric(eventdef(e).baselinewindow)
                    eventdef(e).baselinewindow = eventdef(e).baselinewindow + eventdef(e).trlshift; % in milliseconds
                end
            else
                aas_log(aap,false,['WARNING: timewindow and baselinewindow are exected to be consistent across events\n\t'...
                    '-> only those of the first event will be considered']);
            end
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
            if ~isempty(eventdef)
                epochEEG = pop_epoch(segEEG,{eventdef.eventvalue},eventdef(1).eventwindow);
            else
                epochEEG = segEEG;
            end
                        
            % baseline correction
            if ~isempty(eventdef(1).baselinewindow)
                if isnumeric(eventdef(1).baselinewindow)
                    epochEEG = pop_rmbase(epochEEG,eventdef(1).baselinewindow);
                else
                    switch eventdef(1).baselinewindow
                        case 'all'
                            epochEEG = pop_rmbase(epochEEG,[]);
                        otherwise
                            % NYI
                    end
                end
            end

            datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),sprintf('epochs%s_seg-%d.set',condfn,seg));
            datafn{end+1} = fullfile(aas_getsesspath(aap,subj,sess),sprintf('epochs%s_seg-%d.fdt',condfn,seg));
            pop_saveset(epochEEG,datafn{end});
        end
        
        EL.unload;
        %% Describe outputs
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
    aas_log(aap,false,'Eventvalue cannot be parsed')
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