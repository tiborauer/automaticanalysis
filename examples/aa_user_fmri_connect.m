% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This is an example how to connect fmri to a sructural pipeline to obtain DARTEL normalisation parameters
% It connects to a fully preprocessed structural dataset in /imaging/xy00/aa/Structural
% It uses a tasklist aap_tasklist_fmri_connect.xml (included)
%
% Further features resting state analysis:
%	- Cleaning spikes with BrainWavelet
%	- Cleaning tissue specific global signals
% 	- Temporal filtering
% 	- RSFA estimation (employing aamod_maths)
%
% Tibor Auer, MRC-CBSU
% 01-02-2016

%% INITIALISE
clear

SUBJ = {...
     'S01' 140905; ...
     'S02' 140910; ...
     'S03' 140913; ...
     'S04' 140928; ...
     'S05' 140931; ...
     };
 
aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_fmri_connect.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.NIFTI4D = 1;
aap.options.email='xy00@mrc-cbu.cam.ac.uk';
aap.tasksettings.aamod_slicetiming.sliceorder=[32:-1:1];
aap.tasksettings.aamod_slicetiming.refslice = 16; 
aap.tasksettings.aamod_norm_write_dartel.vox = [3 3 3];
aap.tasksettings.aamod_norm_write_meanepi_dartel.vox = [3 3 3];
aap.tasksettings.aamod_waveletdespike.maskingthreshold=0.7;
aap = aas_renamestream(aap,'aamod_mask_fromsegment_00001','reference','meanepi');
aap = aas_renamestream(aap,'aamod_mask_fromsegment_00001','grey','normalised_grey');
aap = aas_renamestream(aap,'aamod_mask_fromsegment_00001','white','normalised_white');
aap = aas_renamestream(aap,'aamod_mask_fromsegment_00001','csf','normalised_csf');
for b = 1:3
    aap.tasksettings.aamod_firstlevel_model(b).allowemptymodel=1;
    aap.tasksettings.aamod_firstlevel_model(b).writeresiduals=NaN;
    aap.tasksettings.aamod_firstlevel_model(b).moveMat=[1 1 0;1 1 0];

    aap.tasksettings.aamod_temporalfilter(b).filter.LowCutoffFreqency=0.009;
    aap.tasksettings.aamod_temporalfilter(b).filter.HighCutoffFreqency=0.1;
    
    aap.tasksettings.aamod_maths(b).operation='std(X,[],4)';

    aap=aas_renamestream(aap,sprintf('aamod_maths_%05d',b),'input','epi');
    aap=aas_renamestream(aap,sprintf('aamod_maths_%05d',b),'output','RSFA','output');
end

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/aa'; 
aap.directory_conventions.analysisid = 'Functional'; 
aap.directory_conventions.subject_directory_format = 3;

% Add data
aap = aas_addsession(aap,'resting');
aap = aas_addsession(aap,'avtask');
aap = aas_addsession(aap,'movie');
for subj = 1:size(SUBJ,1)
    ser_epi = basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{subj,2},1),'dir','.*CBU_EPI_restingstate$')); ser_epi = deblank(ser_epi(end,:));
    ser{1} = sscanf(ser_epi,aap.directory_conventions.seriesoutputformat);
    ser_epi = basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{subj,2},1),'dir','.*CBU_EPI_sensorimotor_task$')); ser_epi = deblank(ser_epi(end,:));
    ser{2} = sscanf(ser_epi,aap.directory_conventions.seriesoutputformat);
    ser_epi = basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{subj,2},1),'dir','.*CBU_MEPI5_movie$'));
    if isempty(ser_epi) % Prisma
        ser_epi = basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{subj,2},1),'dir','.*Mag_CombinedCoils$')); ser_epi = ser_epi(end-4:end,:);
    end
    for sess = 1:size(ser_epi,1)
        ser{3}(sess) = sscanf(deblank(ser_epi(sess,:)),aap.directory_conventions.seriesoutputformat);
    end
    aap = aas_addsubject(aap,SUBJ{subj,1},SUBJ{subj,2},'functional',ser);
end

%% CONNECT
aap=aas_doprocessing_initialisationmodules(aap);
aap.directory_conventions.allowremotecache = 0;
remotePipes = struct('host',           '', ...
    'directory',      fullfile(aap.acq_details.root,'Structural'), ...
    'allowcache',     0, ...
    'maxstagetag',   'aamod_dartel_normmni_00001', ...
    'checkMD5',       1);
aap=aas_connectAApipelines(aap,remotePipes);

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));