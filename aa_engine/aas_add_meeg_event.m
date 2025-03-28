function aap=aas_add_meeg_event(aap,modulename,subject,session,eventname,eventdef,trialshift,eventwindow,baselinewindow)
% Op.1a.: Using onsets from data --> may multiple calls per subject/session
%     argument 5 (eventname)      - conditionlabel
%     argument 6 (eventvalue)     - eventvalue
%     argument 7 (trialshift)     - trlshift (must not be empty!!!)
%     argument 8 (eventwindow)    - eventwindow ([1x2] array; in ms; must not be empty!!!)
%     argument 9 (baselinewindow) - eventwindow ([1x2] array; in ms; can be empty)
% Op.1b.: Using onsets from data with eventtype specified --> may multiple calls per subject/session
%     argument 5 (eventname)                    - conditionlabel
%     argument 6 (evenettype and eventvalue)    - {evenettype eventvalue}
%     argument 7 (trialshift)                   - trlshift (must not be empty!!!)
%     argument 8 (eventwindow)                  - eventwindow ([1x2] array; in ms; must not be empty!!!)
%     argument 9 (baselinewindow)               - eventwindow ([1x2] array; in ms; can be empty)
% Op.2.: Using user-specific onsets --> single call per subject/session
%     argument 5 (eventname)      - conditionlabels or empty if provided in mat-file along with trl
%     argument 6 (eventvalue)     - trl matrix or mat-file containing trl (and conditionlabels)
%     argument 7 (trialshift)     - []
%     argument 8 (eventwindow)    - eventwindow ([1x2] array; in ms; must not be empty!!!)
%     argument 9 (baselinewindow) - eventwindow ([1x2] array; in ms; can be empty)
% Op.3.: Segment data using onsets from data --> may multiple calls per subject/session
%     argument 5 (eventname)    - 'segment-<segmentlabel>'
%     argument 6 (eventvalue)   - '<eventvalue>:<eventvalue>' or [1x2] array of sample numbers
%     argument 7 (trialshift)   - trlshift (must not be empty!!!)
% If eventtype is empty (also in the module settings!), we assume there is only one eventtype

% Regexp for number at the end of a module name, if present in format _%05d (e.g, _00001)
m1 = regexp(modulename, '_\d{5,5}$');

% Or, we could use '_*' at the end of the module name to specify all modules with that name
m2 = regexp(modulename, '_\*$');

% Or, we might specify specific modules with  '_X/X/X'
m3 = regexp(modulename, '[_/](\d+)', 'tokens');

if ~isempty(m1)
    moduleindex = str2num(modulename(m1+1:end));
    modulename = modulename(1:m1-1);
    
elseif ~isempty(m2)
    modulename = modulename(1:m2-1);
    moduleindex = 1:length(find(strcmp({aap.tasklist.main.module.name}, modulename)));
    
elseif ~isempty(m3)
    modulename = modulename(1:find(modulename=='_',1,'last')-1);
    moduleindex = cellfun(@str2num, [m3{:}]);
    
else
    moduleindex = 1;
end
event = struct('conditionlabel',eventname,...
               'trlshift',trialshift,...
               'eventtype', [],...
               'eventwindow',[],...
               'baselinewindow',[]);
if iscell(eventdef) % specify eventtype
    event.eventtype = eventdef{1};
    event.eventvalue = eventdef{2};
else
    event.eventvalue = eventdef;
end
if ~startsWith(event.conditionlabel,'segment')
    event.eventwindow = eventwindow;
    event.baselinewindow = baselinewindow;
end

% segment definition
if find(strcmp(strsplit(event.conditionlabel,'-'),'segment'),1,'first') == 1
    if numel(event.trlshift) < 2
        aas_log(aap,false,'Segment definition requires trialshift specified as 1x2 array of non-zero values to avoid overlap -> default [2 -2] will be applied');
        event.trlshift = [2 -2];
    end    
end

% find models that corresponds and add events if they exist
for m = 1 : length(moduleindex)
    
    mInd = moduleindex(m);
    
    whichmodel=[strcmp({aap.tasksettings.(modulename)(mInd).condition.subject},subject)] & [strcmp({aap.tasksettings.(modulename)(mInd).condition.session},session)];
    if (~any(whichmodel))
        emptymod=aap.tasksettings.(modulename)(mInd).condition(1);
        emptymod.subject=subject;
        emptymod.session=session;
        if ~isempty(event.trlshift) ... % data-specified
                && isempty(event.eventtype) 
            event.eventtype = aap.tasksettings.(modulename)(mInd).eventtype;
        end
        emptymod.event = event;
        aap.tasksettings.(modulename)(mInd).condition(end+1)=emptymod;
    else
        if ~isempty(event.trlshift) % data-specified
            if isempty(event.eventtype), event.eventtype = aap.tasksettings.(modulename)(mInd).eventtype; end
            aap.tasksettings.(modulename)(mInd).condition(whichmodel).event(end+1) = event;
        else % user-specified
            aas_log(aap,false,sprintf('WARNING: event(s) already specified for %s, %s...Overwriting!!!',subject, session));
            aap.tasksettings.(modulename)(mInd).condition(whichmodel).event = event;
        end
    end
    
end
end