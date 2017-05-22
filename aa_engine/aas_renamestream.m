function aap=aas_renamestream(aap,stage,origstream,newstream,type)
% AAS_RENAMESTREAM 
% Rename a TYPE (input|output) of stream of a STAGE from ORIGSTREAM to NEWSTREAM
%
% E.g. aap = aas_renamestream(aap,'aamod_firstlevel_model_00001','epi','aamod_coreg_extended_2epi_00001.epi')

if nargin < 5
    type = 'input';
end

if strcmp(type,'output')
    aas_log(aap,false,sprintf('WARNING: Renaming output %s of %s to %s.\n  Make sure that you know what you are doing!',origstream,stage,newstream));
end

%% Locate STAGE
[stagename, index] = strtok_ptrn(stage,'_0');
stageindex = sscanf(index,'_%05d');
if ~isfield(aap.tasksettings,stagename) || (numel(aap.tasksettings.(stagename)) < stageindex)
    aas_log(aap,1,sprintf('ERROR: Stage %s not found!',stage));
end
modulenumber = cell_index({aap.tasklist.main.module.name},stagename); modulenumber = modulenumber(stageindex);

%% Check task
odot = origstream == '.';
if any(odot), origstreamname = origstream(find(odot)+1:end);
else origstreamname = origstream; end

ndot = newstream == '.';
if any(ndot), newstreamname = newstream(find(ndot)+1:end);
else newstreamname = newstream; end

toSpecify = any(ndot);
toRename = ~strcmp(origstreamname,newstreamname);

%% Locate ORIGSTREAM
inputstreams = aap.tasksettings.(stagename)(end).([type 'streams']).stream;
if ~iscell(inputstreams)
    inputstreams = {inputstreams};
end
for i = 1:numel(inputstreams)
    idot = inputstreams{i} == '.';
    if any(idot), inputstreamname = inputstreams{i}(find(idot)+1:end);
    else inputstreamname = inputstreams{i}; end
    if strcmp(inputstreamname,origstreamname)
        break;
    end
end
if strcmp(inputstreamname,origstreamname)
    ind = i;
else
    aas_log(aap,1,sprintf('ERROR: Stream %s of module %s not found!',origstream,stagename));
end
stream = aap.schema.tasksettings.(stagename)(end).([type 'streams']).stream{ind};
switch type
    case 'input'
        internalstore = 'inputstreamsources';
    case 'output'
        internalstore = 'outputstreamdestinations';
end
if isfield(aap,'internal')
    storedstreams = aap.internal.(internalstore){modulenumber}.stream;
    if isempty(storedstreams), streamindex = NaN;
    else streamindex = cell_index({storedstreams.name},origstreamname); end
end
        
%% Check ORIGSTREAM
if (~isstruct(stream) || ~isfield(stream.ATTRIBUTE,'isrenameable') || ~stream.ATTRIBUTE.isrenameable) && toRename % no attributes --> assume non-renameable --> allow specifying only
    aas_log(aap,1,sprintf('ERROR: Stream %s of module %s does not support renaming!\nERROR: Check attribute "isrenameable"!',origstream,stagename));
end
if toSpecify
    sourcestage = newstream(1:find(ndot)-1);
    [sourcestagename, index] = strtok_ptrn(sourcestage,'_0');
    sourcestageindex = sscanf(index,'_%05d');
    if ~isfield(aap.tasksettings,sourcestagename) || (numel(aap.tasksettings.(sourcestagename)) < sourcestageindex)
        aas_log(aap,0,sprintf('WARNING: Stage %s not found in local pipeline!',sourcestage));
    end
end

%% Rename ORIGSTREAM
if ~iscell(aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream)
    aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream = newstream;
    aap.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream = newstream;
    if isfield(aap,'internal')
        aap.internal.aap_initial.tasksettings.(stagename)(stageindex).([type 'streams']).stream = newstream;
        aap.internal.aap_initial.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream = newstream;
    end
else
    aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    aap.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    if isfield(aap,'internal')
        aap.internal.aap_initial.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
        aap.internal.aap_initial.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    end
end
if isfield(aap,'internal') && ~isnan(streamindex)
    aap.internal.(internalstore){modulenumber}.stream(streamindex).name = newstream;
    aap.internal.aap_initial.internal.(internalstore){modulenumber}.stream(streamindex).name = newstream;
end
if ~isstruct(aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind})
    aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    if isfield(aap,'internal'), aap.internal.aap_initial.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream; end
else % with attributes
    aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind}.CONTENT = newstream;
    if isfield(aap,'internal'), aap.internal.aap_initial.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind}.CONTENT = newstream; end
end
