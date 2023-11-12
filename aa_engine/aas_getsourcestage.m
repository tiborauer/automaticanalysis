function k = aas_getsourcestage(aap,sourcemod,stream)

currmod = '';
currmodnum = aap.tasklist.currenttask.modulenumber;

if nargin == 2 % based on branch
    currbranch = aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix;
    for currmodnum = find(strcmp({aap.tasklist.main.module.name}, sourcemod))
        srcbranch = aap.tasklist.main.module(currmodnum).extraparameters.aap.directory_conventions.analysisid_suffix;
        if strcmp(currbranch,srcbranch) || contains(currbranch,[srcbranch '_']), break; end % might be an issue if one branchname contains the other
    end
    if ~strcmp(currbranch,srcbranch) && ~contains(currbranch,[srcbranch '_']), currmodnum = aap.tasklist.currenttask.modulenumber; end
    
else % based on stream
    stream = strsplit(stream,'.');
    if numel(stream) > 1
        currmodnum = find(strcmp(arrayfun(@(s) sprintf('%s_%05d',s.name, s.index), aap.tasklist.main.module, 'UniformOutput',false),stream{1}));
        if ~any(strcmp({aap.internal.outputstreamdestinations{currmodnum}.stream.name},stream{2}))
            aas_log(aap,true,sprintf('Stage "%s" does not create stream "%s"',stream{1},stream{2}));
        end
    else
        stream = stream{1};
        while ~strcmp(currmod,sourcemod)
            inpstreams = aap.internal.inputstreamsources{currmodnum}.stream;
            srcmatch = cellfun(@(ss) any(strcmp(ss,stream)), cellfun(@(s) strsplit(s,'.'), {inpstreams.name}, 'UniformOutput', false)); % handle fully specified input streams
            if any(srcmatch)
                currmodnum = inpstreams(srcmatch).sourcenumber;
            else
                aas_log(aap,true,sprintf('stream "%s" cannot be tracked beyond module "%s"',stream,currmod)); 
                currmodnum = [];
                break;
            end
            currmod = aap.tasklist.main.module(currmodnum).name;
        end
    end
end

k = currmodnum;
