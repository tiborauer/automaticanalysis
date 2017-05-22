function pool = qsubcellfun(varargin)
%% Parse input
func = varargin{1};
varargin{1} = func2str(func);
ind_mem = cell_index(varargin,'memreq');
ind_time = cell_index(varargin,'timreq');
ind_stack = cell_index(varargin,'stack');
ind_par = [ind_mem ind_time ind_stack]; ind_par = ind_par(ind_par~=0);
if isempty(ind_par)
    ind_args = 2:nargin;
else
    ind_args = 2:min(ind_par)-1;
end

% Check wheter it is an aa joblist
isaa = isstruct(varargin{2}{1}) && isfield(varargin{2}{1},'options') && ...
    isstruct(varargin{2}{1}.options) && isfield(varargin{2}{1}.options,'aa_minver');
if isaa, aap = varargin{2}{1};
else
    par_xml = spm_select(1,'xml','Select your local parameters and defaults xml...',...
        {fullfile(fileparts(mfilename('fullpath')),'..','aa_recipes_and_parametersets','aap_parameters_defaults_CBSU.xml')},...
        fullfile(fileparts(mfilename('fullpath')),'..','aa_recipes_and_parametersets'));
    aap = xml_read(par_xml,struct('ReadAttr',0));
end

%% Initialise engine
if isaa
    qsubpath = fullfile(getenv('HOME'),'aaworker');
else
    if isempty(aap.directory_conventions.poolprofile)
        aas_log(aap,true,sprintf('poolprofile is not specified in %s',par_xml))
    end
    qsubpath = pwd;
end
qsubpath = [qsubpath filesep func2str(func) '_' datestr(now,30)];
aas_makedir(aap,qsubpath);

try
    profiles = parallel.clusterProfiles;
    if ~any(strcmp(profiles,aap.directory_conventions.poolprofile))
        ppfname = which(spm_file(aap.directory_conventions.poolprofile,'ext','.settings'));
        if isempty(ppfname)
            aas_log(aap,true,sprintf('ERROR: settings for pool profile %s not found!',aap.directory_conventions.poolprofile));
        else
            pool=parcluster(parallel.importProfile(ppfname));
        end
    else
        aas_log(aap,false,sprintf('INFO: pool profile %s found',aap.directory_conventions.poolprofile));
        pool=parcluster(aap.directory_conventions.poolprofile);
    end
    switch class(pool)
        case 'parallel.cluster.Torque'
            aas_log(aap,false,'INFO: Torque engine is detected');
            pool.ResourceTemplate = sprintf('-l nodes=^N^,mem=%dGB,walltime=%d:00:00', aap.options.aaparallel.memory,aap.options.aaparallel.walltime);
        case 'parallel.cluster.Generic'
            aas_log(aap,false,'INFO: Generic engine is detected');
            pool.IndependentSubmitFcn = SetArg(pool.IndependentSubmitFcn,'walltime',aap.options.aaparallel.walltime);
            pool.IndependentSubmitFcn = SetArg(pool.IndependentSubmitFcn,'memory',aap.options.aaparallel.memory);
    end
    pool.NumWorkers = aap.options.aaparallel.numberofworkers;
    pool.JobStorageLocation = qsubpath;
catch ME
    aas_log(aap,false,'ERROR: Cluster computing is not supported!');
    aas_log(aap,true,sprintf('ERROR in %s:\n  line %d: %s',ME.stack.file, ME.stack.line, ME.message));
end

%% Make workers self-sufficient by passing them the paths.
% Users don't need to remember to update
% their own default paths (e.g. for a new aa version)
if isaa
    % AA
    mfp=textscan(which('aaq_qsub'),'%s','delimiter',filesep); mfp = mfp{1};
    mfpi=find(strcmp('aa_engine',mfp));
    aapath=textscan(genpath([filesep fullfile(mfp{1:mfpi-1})]),'%s','delimiter',':'); aapath = aapath{1};
    % SPM
    aapath{end+1}=fileparts(which('spm')); % SPM dir
    p = textscan(path,'%s','delimiter',':'); p = p{1};
    p_ind = cell_index(p,aapath{end}); % SPM-related dir
    for ip = p_ind
        aapath{end+1} = p{ip};
    end
    if isfield(aap.directory_conventions,'spmtoolsdir') && ~isempty(aap.directory_conventions.spmtoolsdir)
        SPMTools = textscan(aap.directory_conventions.spmtoolsdir,'%s','delimiter', ':');
        SPMTools = SPMTools{1};
        for p = SPMTools'
            if exist(p{1},'dir'), aapath{end+1}=p{1};end
        end
    end
    % MNE
    if isfield(aap.directory_conventions,'mnedir') && ~isempty(aap.directory_conventions.mnedir)
        if exist(fullfile(aap.directory_conventions.mnedir,'matlab'),'dir')
            aapath{end+1}=fullfile(aap.directory_conventions.mnedir,'matlab','toolbox');
            aapath{end+1}=fullfile(aap.directory_conventions.mnedir,'matlab','examples');
        end
    end
    aapath=aapath(strcmp('',aapath)==0);
else
    aapath ={};
end

%% Submit
for iJob = 1:numel(varargin{2})
    
    pause(0.5); % do not overload
    
    J = createJob(pool);
    inparg = {};
    nArg = 0;
    for iArg = ind_args
        nArg = nArg + 1;
        inparg{nArg} = varargin{iArg}{iJob};
    end
    
    if isprop(J,'AdditionalPaths')
        J.AdditionalPaths = aapath;
    elseif isprop(J,'PathDependencies')
        J.PathDependencies = aapath;
    end
    
    createTask(J,func,0,inparg);
    fprintf('SUBMIT: %s\n',func2str(func));
    
    J.submit;
end

end

function argout = SetArg(argin,key,value)
argout = argin;
if ~iscell(argout), argout = {argout}; end
ind = find(cellfun(@(x) strcmp(x,key),argout));
if ind, argout(ind:ind+1) = []; end
argout(end+1:end+2) = {key value};
end