function aap = aas_addsubject(aap, varargin)
% Add subject/data to the analysis. It can be called multiple times to add more subjects and/or more sources to a particular subject.
%
% FORMAT function aap = aas_addsubject(aap, data)
% Process only autoidentified images. User will be a warned about the lack of series specification. Subject name, which is used as a reference in aa 
% (i.e. when calling aas_addevent, aas_addcontrasts, etc.), will be automatically defined based on aap.directory_conventions.subject_directory_format:
%   0   - based on predefined list stored in aap.directory_conventions.subject_directory_names
%   1   - based on the "data" (see below)
%   2   - based on the order of specification (S01, S02, etc.)
%
% aap           - aap structure with parameters and tasklist
% data          - subject foldername within database. 
%                   - for MRI: a single entry according to aap.directory_conventions.subjectoutputformat
%                   - for MEEG: it is a cell array of two entries according to aap.directory_conventions.meegsubjectoutputformat (1st entry for MEEG data)
%                     and aap.directory_conventions.subjectoutputformat (2nd entry for MRI data). When MRI data is not analysed, the 2nd entry must be 
%                     an empty array.
%
%
% FORMAT function aap = aas_addsubject(aap, name, data)
% As above, but manually specifies subject name. aap.directory_conventions.subject_directory_format must be set to 3.
%
% name          - subject name
%
%
% FORMAT function aap = aas_addsubject(___,'name',subjectname)
% Another way to specify subject name manually, if aap.directory_conventions.subject_directory_format is not set to 3.
%
% subjectname   - subject name as text string
%
%
% FORMAT function aap = aas_addsubject(___,'functional',series)
% Specify functional (fMRI/MEEG) data.
%
% series        - for DICOM: array of series number(s) of EPIs. E.g.: 
%                   two series of single-echo EPI: [5 10]
%                   two series of single-echo EPI and one series of multi-echo EPI with 5 echos: {5 10 15:19}
%               - for NIfTI: cell array containing one or more
%                   for structural: string containing a full or relative path (from a single rawdatadir)
%                   for 4D NIfTI: string containing a full or relative path (from a single rawdatadir)
%                   for whole-brain EPI: string containing a full or relative path (from a single rawdatadir). Can be specified only after fMRI series.
%                   for 3D NIfTI: cell array (i.e. nested) of strings of full path
%                 Strings can be replaced by structures with fields 'fname' (path to image) and 'hdr' (path to header) to specify metadata.
%               - for MEEG: full or relative filename of the acquisition file in the subject folder
% Series have to be specified in the same order as the corresponding sessions have been added in the UMS. Missing series can be specified either with "0" (for numerical array input) or with "[]" (for cell array input).
%
%
% FORMAT function aap = aas_addsubject(___,'diffusion',series)
% Specify diffusion-weighted MRI data.
%
% series        - for DICOM: numeric array of series number(s)
%               - for NIfTI: cell of structure(s) with fields 'fname' (path to image), and 'bval', 'bvec'(path to bvals and bvecs)
% Series have to be specified in the same order as the corresponding sessions have been added in the UMS. 
%
%
% FORMAT function aap = aas_addsubject(___,'structural', series)
% Specify structural data (overwrites autoidentification).
%
% series        - for DICOM: numeric array of series number
%               - for NIfTI: cell containing a string (path to image) or a structure with fields 'fname' (path to image) and 'hdr' (path to header)
%
%
% FORMAT function aap = aas_addsubject(___,'fieldmaps', series)
% Specify fieldmap data (overwrites autoidentification).
%
% series        - for DICOM: numeric array of series numbers
%               - for NIfTI: cell of structure with fields 'fname' (cell of 3 filenames - 2x magnitude + 1x phase), 'hdr' (path to header), and 'session' (cell of session names or '*' for all sessions)
%
%
% FORMAT function aap = aas_addsubject(___,'specialseries', series)
% Specify 'special' data (e.g. ASL, MTI, MPM).
%
% series        - for DICOM: cell array of numeric arrays of series numbers
%               - for NIfTI: not supported yet
% Series have to be specified in the same order as the corresponding sessions have been added in the UMS. 
%
%
% FORMAT function aap = aas_addsubject(___,'ignoreseries', series)
% Specify DICOM series to be ignored during autoidentification.
%
% series        - numeric arrays of series numbers

%% Parse
iMRIData = 1; % new subject
iMEEGData = 1;
if isempty(aap.acq_details.subjects(end).subjname)
    subjind = 1;
else
    subjind = numel(aap.acq_details.subjects) + 1;
end

name = '';
switch aap.directory_conventions.subject_directory_format
    case 0 % from predefined list
        name = aap.directory_conventions.subject_directory_names{subjind};
        data = varargin{1};
        varargin(1)= [];
    case 1 % from data
        data = varargin{1};
        varargin(1)= [];
    case 2 % S#
        name = sprintf('S%02d',subjind);
        data = varargin{1};
        varargin(1)= [];
    case 3 % manual
        name = varargin{1};
        data = varargin{2};
        varargin(1:2)= [];
    otherwise
        aas_log(aap,true,sprintf('ERROR: Unknown subject directory format (aap.directory_conventions.subject_directory_format=%d. Value only 0-3 is allowed.',aap.directory_conventions.subject_directory_format));
end
try
    args = vargParser(varargin);
catch
    aas_log(aap,false,sprintf('ERROR in %s: incorrect arguments',mfilename),'Errors')
    help(mfilename);
    error('ERROR in %s: incorrect arguments',mfilename);
end

%% Sanity, compatiblity check
if isempty(varargin)
    aas_log(aap,false,'WARNING: No series has been specified!\n')
else
    if ~isa(varargin{1},'char')
        aas_log(aap,true,sprintf('ERROR: Arguments are  different from what expected!\n %s',help('aas_addsubject')))
    end
end

%% Initialize subject
% with a blank template for a subject entry
fields=fieldnames(aap.schema.acq_details.subjects);
fields(strcmp(fields,'ATTRIBUTE')) = [];
for field=fields'
    thissubj.(field{1})={[]};
end
fields(strcmp(fields,'subjname')) = [];

% search for existing subject
if isfield(args,'name'), name = args.name; end
if ~isempty(name) && ~isempty(aap.acq_details.subjects(1).subjname)
% name specified --> check whether subject already exists (only if there is at least one already)    
    subjserach = cell_index({aap.acq_details.subjects.subjname},name);
    if subjserach
        subjind = subjserach; 
        thissubj = aap.acq_details.subjects(subjind);
        iMRIData = numel(thissubj.mriname)+1;
        iMEEGData = numel(thissubj.meegname)+1;
        for field=fields'
            thissubj.(field{1}){end+1}=[];
        end
    end
end

%% Data
try
    if iscell(data) && numel(data) == 2 % MEEG
        thissubj.meegname{iMEEGData}=data{1};
        thissubj.mriname{iMRIData}=data{2};
        if isempty(name), name = aas_meegname2subjname(aap,sprintf(aap.directory_conventions.meegsubjectoutputformat,thissubj.meegname{1})); end
    else % MRI
        thissubj.mriname{iMRIData}=data;
        if isempty(name), name = aas_mriname2subjname(aap,sprintf(aap.directory_conventions.subjectoutputformat,thissubj.mriname{1})); end
    end
catch
    aas_log(aap,true,'In aas_addsubject, name is expected to be either single string for MRI, or a cell of two for MEEG written like this {''meegname'',''mriname''}.');
end

thissubj.subjname = name;

%% Series
if isfield(args,'functional') && ~isempty(args.functional) 
    if isnumeric(args.functional) || isnumeric(args.functional{1}) % DICOM series number --> MRI
        thissubj.seriesnumbers{iMRIData}=args.functional;
    else
        fMRI = {}; MEEG = {};
        fMRIdim = [];
        for s = 1:numel(args.functional)
            if iscell(args.functional{s}) % multiple 3D files
                fMRI{end+1} = args.functional{s};
            elseif isempty(args.functional{s}) % missing series
                fMRI{end+1} = [];
                MEEG{end+1} = [];
            elseif ischar(args.functional{s}) ||... % NIfTI file
                    isstruct(args.functional{s})    % hdr+fname
                % Get filename
                
                if isstruct(args.functional{s})
                    if numel(args.functional{s}.fname) > 1 % multiple 3D files
                        fMRI{end+1} = args.functional{s};
                        continue;
                    end
                    fname = args.functional{s}.fname;
                else
                    fname = args.functional{s};
                end
                
                % - try in rawmeegdatadir
                if ~exist(fname,'file')
                    if ~isempty(thissubj.meegname{iMEEGData})
                        tmpaap = aap;
                        tmpaap.directory_conventions.meegsubjectoutputformat = '%s';
                        if exist(fullfile(meeg_findvol(aap,thissubj.meegname{iMEEGData},'fullpath',true),fname),'file') ||...
                                ~isempty(meeg_findvol(tmpaap,fname)) % try empty room
                            MEEG{end+1} = fname;
                            continue;
                        end
                    end
                end
                % - try in rawdatadir
                if ~exist(fname,'file'), fname = fullfile(aas_findvol(aap,spm_file(fname,'path')),spm_file(fname,'filename')); end
                
                if ~exist(fname,'file'), aas_log(aap,1,sprintf('ERROR: File %s does not exist!',fname)); end
                
                % Sort
                if strcmp(spm_file(fname,'ext'),'fif')
                     MEEG{end+1} = fname;
                     continue;
                end
                
                V = spm_vol(fname);
                if numel(V) > 1 % 4D --> fMRI
                    fMRI{end+1} = args.functional{s};
                    fMRIdim = V(1).dim;
                else % 3D --> structural
                    if ~isempty(fMRIdim) && all(fMRIdim(1:2) == V.dim(1:2)) % same inplane resolution as fMRI
                        thissubj.wholebrain_epi{iMRIData}=args.functional(s);
                    else
                        thissubj.structural{iMRIData}=args.functional(s);
                    end
                end
            else % mixed: DICOM series number for fMRI
                thissubj.seriesnumbers{iMRIData}=args.functional{s};
            end
        end
        if ~isempty(fMRI) && any(cellfun(@(x) ~isempty(x), fMRI))
            thissubj.seriesnumbers{iMRIData}=fMRI;
        end
        if ~isempty(MEEG) && any(cellfun(@(x) ~isempty(x), MEEG))
            thissubj.meegseriesnumbers{iMEEGData}=MEEG;
        end
    end
end

if isfield(args,'diffusion') && ~isempty(args.diffusion)
    thissubj.diffusion_seriesnumbers{iMRIData}=args.diffusion;
end

for meas = {'structural' 'fieldmaps' 'specialseries' 'ignoreseries'}
    if isfield(args,meas{1}) && ~isempty(args.(meas{1}))
        thissubj.(meas{1}){iMRIData}=args.(meas{1});
    end
end

% And put into acq_details, replacing a single blank entry if it exists
aap.acq_details.subjects(subjind)=thissubj;
end
