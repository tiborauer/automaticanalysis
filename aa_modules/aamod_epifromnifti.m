% AA module - EPI from NIFTI

function [aap,resp]=aamod_epifromnifti(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        %% Init
        
        numdummies = aap.acq_details.numdummies;
        if ~isempty(aap.tasklist.currenttask.settings.numdummies)
            numdummies=aap.tasklist.currenttask.settings.numdummies;
        end
        
        %% Select
        series = horzcat(aap.acq_details.subjects(subj).seriesnumbers{:});
        if ~iscell(series) ... 
                || (~isstruct(series{1}) ... % hdr+fname
                && ~ischar(series{1}) ... % fname
                && ~iscell(series{1})) % fname (4D)
            aas_log(aap,true,['ERROR: Was expecting list of struct(s) of fname+hdr or fname in cell array\n' help('aas_addsubject')]);            
        end
        series = series{sess};
        
        %% Process
        sesspth=aas_getsesspath(aap,subj,sess);
        aas_makedir(aap,sesspth);

        % Files
        headerFn ='';
        imageFn = series;
        if isstruct(imageFn)
            headerFn = imageFn.hdr;
            imageFn = imageFn.fname;
        end
        headerfile = headerFn;
        niftifile = imageFn;
        if ~iscell(niftifile), niftifile = {niftifile}; end % 4D-NIFTI
        if ~exist(niftifile{1},'file') % assume path realtive to the subject's dir
            niftisearchpth=aas_findvol(aap,subj);
            if ~isempty(niftisearchpth)
                niftifile = cell_fullfile(niftisearchpth,niftifile);
                if ~exist(niftifile{1},'file'), aas_log(aap,true,['ERROR: image ' niftifile{1} ' not found']); end
                if ~isempty(headerFn)
                    headerfile = fullfile(niftisearchpth,headerFn);
                end
            end
        end
        comp = false;
        if any(strcmp(spm_file(niftifile,'Ext'),'gz'))
            comp = true;
            gunzip(niftifile{1},fullfile(sesspth,'temp'));
            niftifile{1} = spm_file(niftifile{1}(1:end-3),'path',fullfile(sesspth,'temp'));
        end
        
        % Header
        DICOMHEADERS{1} = struct;
        if ~isempty(headerfile)
            if strcmp(spm_file(headerfile,'Ext'),'mat')
                dat = load(headerfile);
                DICOMHEADERS = dat.DICOMHEADERS;
            else
                switch spm_file(headerfile,'Ext')
                    case 'dcmhdr'
                        DICOMHEADERS{1} = header_dcmhdr(headerfile);
                    case 'json'
                        DICOMHEADERS{1} = header_json(headerfile); % BIDS
                end
                DICOMHEADERS{1}.volumeTR = DICOMHEADERS{1}.RepetitionTime/1000;
                DICOMHEADERS{1}.volumeTE = DICOMHEADERS{1}.EchoTime/1000;
                DICOMHEADERS{1}.slicetimes = DICOMHEADERS{1}.SliceTiming/1000;
                [junk, DICOMHEADERS{1}.sliceorder] = sort(DICOMHEADERS{1}.slicetimes);
                DICOMHEADERS{1}.echospacing = DICOMHEADERS{1}.EchoSpacing/1000;
            end
        else
            aas_log(aap,false,'WARNING: No header provided!');
        end        
        
        % Image
        finalepis={};
        V = spm_vol(niftifile); 
        if iscell(V), V = cell2mat(V); end
        fle = spm_file(niftifile,'basename');
        ext = spm_file(niftifile,'Ext');
        for fileind=1:numel(V)
            Y=spm_read_vols(V(fileind));
            if numel(niftifile) == 1 % 4D-NIFTI
                fn=fullfile(sesspth,sprintf('%s_%04d.%s',fle{1},fileind,ext{1}));
            else % 3D-NIFTI
                fn=fullfile(sesspth,[fle{fileind} '.' ext{fileind}]);
            end
            V(fileind).fname=fn;
            V(fileind).n=[1 1];
            spm_write_vol(V(fileind),Y);
            finalepis = [finalepis fn];
        end
        
        if isfield(DICOMHEADERS{1},'PhaseEncodingDirection') && ~isempty(DICOMHEADERS{1}.PhaseEncodingDirection)
            sliceaxes = {'x' 'y'};
            ind = cell_index(sliceaxes, DICOMHEADERS{1}.PhaseEncodingDirection(1));
            if ind == 0
                % newer BIDS spec uses this format insteadj
                sliceaxes = {'i' 'j'};
                ind = cell_index(sliceaxes, DICOMHEADERS{1}.PhaseEncodingDirection(1));
                if ind == 0
                    aas_log(aap,1,sprintf(...
                        'Could not parse PhaseEncodingDirection: %s',...
                        DICOMHEADERS{1}.PhaseEncodingDirection(1)));
                end
            end
            DICOMHEADERS{1}.NumberOfPhaseEncodingSteps = V(1).dim(ind);
        end
        
        % Write out the files
        
        % Now move dummy scans to dummy_scans directory
        dummylist=[];
        if numdummies
            dummypath=fullfile(sesspth,'dummy_scans');
            aap=aas_makedir(aap,dummypath);
            for d=1:numdummies
                cmd=['mv ' finalepis{d} ' ' dummypath];
                dummylist=strvcat(dummylist,fullfile(dummypath, spm_file(finalepis{d},'filename')));
                [s, w]=aas_shell(cmd);
                if (s)
                    aas_log(aap,1,sprintf('Problem moving dummy scan\n%s\nto\n%s\n\n%s',finalepis{d},dummypath,w));
                end
            end
        else
            d = 0;
        end
        finalepis = {finalepis{d+1:end}};
        % 4D conversion [TA]
        if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D
            finalepis = finalepis{1};
            ind = regexp(finalepis,'_[0-9]*\.'); % find volume numbers just before the extension
            sfx = '';
            if isempty(ind)
                ind = find(finalepis=='.');
                sfx = '_4D';
            end
            finalepis = [finalepis(1:ind-1) sfx '.nii'];
            if iscell(V), V = cell2mat(V); end
            spm_file_merge(char({V(numdummies+1:end).fname}),finalepis,0,DICOMHEADERS{1}.volumeTR);
        end
              
        % Describe outputs

        % new: output streams are renameable so that other modalities
        % (e.g., dot) can reuse aamod_epifromnifti w/ just a header alias
        
        % we assume the first output streamname defines the modality;
        % the remaining outstream names are derived from it
        
        outputstreams = aas_getstreams(aap,'output'); 
        
        ostream_01 = outputstreams{1};
       
        ostream_02 = 'dummyscans';
        
        if strcmp(ostream_01,'epi')
            ostream_03 = 'epi_dicom_header';
        elseif strcmp(ostream_01,'dot')
            ostream_03 = 'dot_header';
        else
            ostream_03 = 'functional_header';
        end      
         
        if comp, rmdir(fullfile(sesspth,'temp'),'s'); end
        aap=aas_desc_outputs(aap,subj,sess, ostream_01,finalepis);
        aap = aas_desc_outputs(aap,subj,sess,ostream_02,dummylist);
        
        dcmhdrfn = fullfile(sesspth,'functional_header.mat');
        save(dcmhdrfn,'DICOMHEADERS');
        aap = aas_desc_outputs(aap,subj,sess,ostream_03,dcmhdrfn);

        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
end

function hdr = header_dcmhdr(fname)
% minimal fileds
hdr.RepetitionTime = [];
hdr.EchoTime = [];
hdr.PhaseEncodingDirection = '';
hdr.EchoSpacing = [];
hdr.SliceTiming = [];

flines = cellstr(fileread(fname));

inds = strfind(flines,'Repetition Time');
for i = 1:numel(inds)
    if ~isempty(inds{i}), break; end
end
str = textscan(flines{i}(inds{i}:end),'%s');
hdr.RepetitionTime = str2double(str{1}{end});

inds = strfind(flines,'Echo Time');
for i = 1:numel(inds)
    if ~isempty(inds{i}), break; end
end
str = textscan(flines{i}(inds{i}:end),'%s');
hdr.EchoTime = str2double(str{1}{end});
end

function hdr = header_json(fname) % BIDS
% minimal fileds
hdr.RepetitionTime = [];
hdr.EchoTime = [];
hdr.PhaseEncodingDirection = '';
hdr.EchoSpacing = [];
hdr.SliceTiming = [];

% retrieve info
info = loadjson(fname);
for f = fieldnames(info)'
    hdr.(f{1}) = info.(f{1});
end

if isfield(hdr,'EffectiveEchoSpacing')
    hdr.EchoSpacing = hdr.EffectiveEchoSpacing*1000;
end

% convert timings to ms (DICOM default)
hdr.RepetitionTime = hdr.RepetitionTime*1000;
hdr.EchoTime = hdr.EchoTime*1000;
hdr.SliceTiming = hdr.SliceTiming*1000;
end

function fullfilepaths=cell_fullfile(basepath,filepaths)
% Like full file, except that input filepaths is a cell array. basepath is
% prepended (using fullfile) to every item in filepaths. Result is a cell
% array
for fileind=1:length(filepaths)
    fullfilepaths{fileind}=fullfile(basepath,filepaths{fileind});
end
end
