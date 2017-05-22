## 5.3.0 ##

#### New general features ####
  - [aas_processBIDS](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_processBIDS.m) generates subject names without "subj-"
  - [aas_processBIDS](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_processBIDS.m) and [aa_export_toBIDS](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_toolbox/aa_export_toBIDS.m) supports T2-weighted and diffusion data
  - [aaq_qsub](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_qsub.m) and [aaq_matlab_pct](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_matlab_pct.m) are based on the same setting `aap.directory_conventions.poolprofile` which refers to a more standard profile rather than to a custom function
  - module domain and modality is changed dynamically according to the input
  - processing ASL data
  - `ignoreseries` is implemented
  - improved compatibility with MATLAB R2015a
  - improved compatibility with SPM8
  - stand-alone:
    - pragmas added to support building stand-alone version
    - [aa_build_standalone](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_toolbox/aa_build_standalone.m) to generate new functions with actualised pragmas
  - [aamod_math](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_math.xml) supports
    - series of operations
    - arithmetic operations (+,-,*,/)
  
  #### New features for structural ####
  - [aamod_convert_structural](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_convert_structural.xml) supports averaging
  - [aamod_structuralfromnifti](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_structuralfromnifti.xml)
    - handles multiple images (+averaging)
	- handles multiple modalities (T1, T2)
	- saves header (BIDS)
  
  #### New features for Diffusion ####
  - [aamod_convert_diffusion](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_convert_diffusion.xml) outputs DICOM header
  - [aamod_diffusion_coreg_extended_2DKI](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_diffusion_coreg_extended_2DKI.xml) and [aamod_diffusion_coreg_extended_2DTI](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_diffusion_coreg_extended_2DTI.xml) are parallelised across sessions.
  
  #### New features for fMRI ####
  - [aas_processinput](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_processinput.m) can process pipeline without modelling and pipeline with multiple models
  - sessions can be contrasted on first-level
  - [aamod_coreg_extended_2epi](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_coreg_extended_2epi.xml) is parallelised across sessions. N.B.: It will not output meanepi (see new module `aamod_coreg_extended_2meanepi`)
  - new modules:
    - [aamod_coreg_extended_2meanepi](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_coreg_extended_2meanepi.xml)
    - [aamod_split_epi](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_split_epi.xml) for splitting (subset of) a 4D volume (into 3D)

  #### Bugfixes ####
  - [aaq_matlab_pct](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_matlab_pct.m) can launch arbitrary number of workers 
  - handling multiple structurals (e.g. ME-FLASH)
  - [aamod_epifromnifti](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_epifromnifti.xml) handles
    - 3D NIfTI source
	- path relative to `rawdatadir`
	- manual header
  - [aamod_secondlevel_threshold](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_threshold.xml):
	- default threshold is set to p = 0.05 FWE-corrected
    - output for each first-level contrasts
  - [aa_export_toBIDS](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_toolbox/aa_export_toBIDS.m) handles
    - Diffusion data without DICOM header
	- fMRI data without events
  
## 5.2.0 ##

#### New general features ####
  - MTI and sepcial_sessions) has benn moved to session domain and thus can be missing
  - `aamod_get_dicom_*` modules can now also edit DICOM fields
  - bounding box can be specified for `aamod_norm_write*` modules
  - `aamod_norm_write*` modules can process deformations
  - `aamod_roi_extract*` has option to use either the space of the data or the ROI
  - subjects with different sessions can be analysed
    - [aamod_firstlevel_contrasts](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_firstlevel_contrasts.xml) can automatically assign contrasts, unless `uniquebysession` is used. :warning: **Contrast names must be unique both within- and across-sessions!**
  - [aamod_coreg_general](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_coreg_general.m) is prepared for session domain and 4D NIfTI
  - killed and finished jobs will not be deleted 
  - queue viewer is optional
  - aaparallel parameters has been added
  - [aaq_matlab_pct](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_matlab_pct.m) has been tested on MATLAB 2015a
  - [aaq_matlab_pct](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_matlab_pct.m) accepts
    - aaparallel parameters
    - external poolprofile function specified in `aap.directory_conventions.poolprofile`
    - pass `aaworker`
  - TFCE ([external](http://dbm.neuro.uni-jena.de/tfce/)) has been updated

#### New features for fMRI ####
  - NIfTI input without header is now allowed
  - examples [aa_user_branching_example1](https://github.com/rhodricusack/automaticanalysis/blob/v5-devel/examples/aa_user_branching_example1.m) and [aa_user_branching_example2](https://github.com/rhodricusack/automaticanalysis/blob/v5-devel/examples/aa_user_branching_example2.m) have been updated
  - estimate and reslice flags has been added to [aamod_realign](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_realign.xml)
  - new modules:
    - [aamod_reorientto_epi](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_reorientto_epi.xml) and [aamod_reorientto_fieldmap](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_reorientto_fieldmap.xml) for longitudinal studies 
  - tsdiffana ([external](https://sourceforge.net/p/spmtools/tsdiffana/code/ci/master/tree)) has been updated
  - Motion FingerPrint ([external](https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software/)) has been updated
	
#### New features for MEG ####
  - `transform` setting of [aamod_meg_maxfilt](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_meg_maxfilt.xml) can be negative to specify session based on `selected_session`
  - [aamod_meg_average_merged](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_meg_average_merged.xml) can deal with merged data
  - new modules:
    - [aamod_meg_merge](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_meg_merge.xml): merging MEG sessions
    - [aamod_meg_grandmean](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_meg_grandmean.xml): estimating grabd mean across subjects

#### Bugfixes ####
  - [aamod_firstlevel_threshold](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_firstlevel_threshold.xml) and [aamod_secondlevel_model](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_model.xml) now can deal with path with "spm" in it
  - [aamod_firstlevel_threshold](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_firstlevel_threshold.xml) and [aamod_secondlevel_model](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_model.xml) now can deal with cases when reslicing activation maps "vanishes" activations
  - modality detection for MTI-related modules
  - [aamod_autoidentifyseries_siemens](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_autoidentifyseries_siemens.m) with missing series
  - [aamod_autoidentifyseries_siemens](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_autoidentifyseries_siemens.m) with custom `dicom_protocol_field`
  - [aas_getfiles_bystream](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_getfiles_bystream.m) now also uses `aap.options.maximumretry` when copying streams potentially used jointly with other jobs
  - coregister after `aamod_norm_noss*`
  - reporting connected pipeline
  - reporting branched pipeline with `selected_sessions`
  - reporting `aamod_*dartel_denorm*`
  - reporting MTI
  - [aamod_bet](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_bet.xml) and derived modules
  - [aamod_MTI_roi_valid](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_MTI_roi_valid.xml)
  
## 5.1.0 ##

#### New general features ####
  - [aa_export_toBIDS](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_toolbox/aa_export_toBIDS.m) export raw data from aa pipeline in BIDS format (limited functionality)

#### New features for structural ####
  - automatic inputstream detection (and renaming) for aamod_roi_valid_structural

#### New features for fMRI ####
  - automatic inputstream detection (and renaming) for aamod_roi_valid_epi
- full set of secondlevel modules: [aamod_secondlevel_model](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_model.xml), [aamod_secondlevel_contrasts](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_contrasts.xml), [aamod_secondlevel_threshold](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_threshold.xml), [aamod_secondlevel_threshold_register2FS](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_threshold_register2FS.xml)

#### New features for MEG ####
  - downsampling added to [aamod_meg_maxfilt](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_meg_maxfilt.xml)

#### Bugfixes ####
  - economise code: aamod_roi_extract and aamod_roi_valid have been replaced/expanded with [aamod_roi_extract_epi](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_roi_extract_epi.xml), [aamod_roi_extract_structural](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_roi_extract_structural.xml), [aamod_roi_valid_epi](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_roi_valid_epi.xml) and [aamod_roi_valid_structural](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_roi_valid_structural.xml)
  - aas_log works with empty `aap`
  - [aas_runfslcommand](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_runfslcommand.m) and [aas_runFScommand](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_runFScommand.m) now pass run-time MATLAB path to in-shell MATLAB (if applicable)
  - [aa_provenance](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_toolbox/provenance/aa_provenance.m) now uses run-time dependency (see "New features ...")
  - [QueueViewer](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/QueueViewerClass.m) will not delete jobs finished with errors --> job folder will be kept for debugging
  - Termination of the pipeline via [aaq_qsubVeiwerClass](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_qsubVeiwerClass.m) is correctly detected by [aaq_qsub](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aaq/aaq_qsub.m) and [aa_doprocessing](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aa_doprocessing.m)
  - [aamod_fieldmapfromnifti](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_fieldmapfromnifti.xml) now correctly outputs TEs from specified header
  - Reporting includes distributions of (both the first- and second-level) contrasts
  - `chainsearch` and `threshold` options have been implemented in [aamod_waveletdespike](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_waveletdespike.xml)
  - [aamod_firstlevel_threshold](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_firstlevel_threshold.xml) and [aamod_secondlevel_model](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_model.xml) now correctly outputs coronal sections
  - [aamod_firstlevel_threshold](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_firstlevel_threshold.xml) and [aamod_secondlevel_model](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_modules/aamod_secondlevel_model.xml) saves outputs correctly even if path contains "spm"
  - [aas_add_meg_session](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel/aa_engine/aas_add_meg_session.m) now prevents adding the same session multiple times (e.g. when the function is called per subject)
  
## 5.0.0 ([branch v5 initial](https://github.com/rhodricusack/automaticanalysis/tree/v5-devel)) ##

As the change in the major versioning implies, older user master scripts are not compatible with v5. Examples has been updated to demonstrate new syntax. In addition, pipelines processed with older versions cannot be re-processed with v5 (`aap` structure stored in _aap\_parameters.mat_ is not compatible). A script `aa_convert_subjects` is provided to convert `aap` structure stored in _aap\_parameters.mat_. When connecting to a remote pipeline processed with an older version of aa, `aa_convert_subjects` is automatically called; so no explicit conversion is required. 

#### New general features ####
  - More detailed documentation of several key functions
  - Explicit subject identifier `aap.acq_details.subjects.subjname`: `aas_addsubject`
    - Subject name is more unambiguously specified
    - Subject name is not tied to the data
    - Same subject name can be used as a reference in the whole UMS

  - Longitudinal/multi-visit measurement
  - In case of selected_session, only relevant inputs will be retrieved
  - Session-specific fieldmaps
  - `aas_addsubject` has a more intuitive parameterisation
  - From NIfTI
    - [BIDS](bids.neuroimaging.io) datasets are supported (full)
  - Reporting remote pipeline
  - Lightweighting: remove ANTs, VBM8 and FreeSurfer deface templates from the package and mark them as (optional) requirements

#### New features for structural ####
  - aamod_mask_fromsegment accepts different exact thresholds for GM, WM, CSF

#### New features for fMRI ####
  - Reorienting input images (structural, diffusion and EPI) to their middle voxel (`aamod_reorienttomiddle_*`)
  - Specifying contrast for certain sessions using format "sessions:<session name>[+<session name>[...]]"
  - Specifying contrast with condition names in a format <weight>x<regressor name>[<main ('m') or parametric ('p')><number of basis/parametric function>] (e.g. '+1xTASK|-1xREST' or '+1xTASKp1|-1xRESTm1'). N.B.: It requires regressor names with UPPERCASE letters only!
	
#### Bugfixes ####
  - Scaling automatic temporal modulation
  - aamod_bet_meanepi
  - aamod_waveletdespike using explicit brainmask

## 4.3.1 ##

#### Bugfix ####
  - Exmaple user master script and tasklist are provided for [demo dataset](http://cusacklab.org/downloads/aa_demo_v1.tar.gz):  ([aa\_user\_demo\_v2](https://github.com/rhodricusack/automaticanalysis/tree/devel-share/examples/aa_user_demo_v2.m))
  - Compatibility with MATLAB r2015b

## 4.3.0 ([branch devel-share](https://github.com/rhodricusack/automaticanalysis/tree/v4-devel)) ##

#### Copyrights ####
  - Sorting _external_/_extrafunctions_
  - Acknowledging externals in [README.json](https://github.com/rhodricusack/automaticanalysis/tree/v4-devel/external/README.json)

#### New general features ####
  - Dynamic modification of streams (`aas_renamestream`)
  - Selected_sessions now works also for branching

  - From NIfTI
    - 3D NIfTI inputs are fully supported
    - NIfTI fieldmaps are supported (`aamod_fieldmapfromnifti`)
    - [BIDS](bids.neuroimaging.io) datasets are supported (partial)
    - `aamod_fsl_reorienttoMNI`
	
#### New features for fMRI ####
  - New fMRI example to demonstrate some (new) features: ([aa\_user\_fmri\_advanced](https://github.com/rhodricusack/automaticanalysis/tree/v4-devel/examples/aa_user_fmri_advanced.m))
    - Specifying structural series
    - Explaining Motion FingerPrint instead of simple mocoparameters in design
    - Automatic slicetiming with exact timing from DICOM header (`autodetectSO`)
    - DARTEL normalisation
    - Activation maps projected to surface using FreeSurfer (`aamod_firstlevel_threshold_register2FS`)
    - Second-level GIFT (`aamod_secondlevel_GIFT`)

  - automatic temporal modulation
  - `aamod_firstlevel_model`: option to save residual 
  - Despiking with BrainWavelet (`aamod_waveletdespike`)
  - `aamod_maths`
  - `aamod_temporalfilter`

#### New features for MEG ####
  - Diagnostics (courtesy to [Rik Henson](https://www.mrc-cbu.cam.ac.uk/people/rik.henson))
  - Update examples: [aa\_user\_meg.m](https://github.com/rhodricusack/automaticanalysis/tree/v4-devel/examples/aa_user_meg.m), [aa\_user_meg\_connect](https://github.com/rhodricusack/automaticanalysis/tree/v4-devel/examples/aa_user_meg_connect.m)
  - Maxfilter (`aamod_meg_maxfilt`) allows getting HPI from a session
  - Maxfilter (`aamod_meg_maxfilt`) accepts custom calibration files
  - Better delegeation of denoising tasks:
    - `aamod_meg_denoise_ICA_1`: runs ICA
    - `aamod_meg_denoise_ICA_1`: thresholds and removes
  - Epoching: `aas_add_meg_event` and `aamod_meg_epochs`
  - Averaging: `aamod_meg_average`

#### Improving robustness ####
  - Code development
    - Simplifying code
      - Setting modality --> `aas_getfile_bystream` automatically detects sessions based on modality
      - Avoid `eval`
      - Avoid `spm_jobman`
    - Removing 'orphan' functions
    - Update: spm_mods based on SPM12 r6470
				
  - Running
    - Compatibility with MATLAB pre-r2012b (local execution only)
    - `qsub` cleans jobs from previous execution
    - Iterative file retrieval (`aap.options.maximumretry`) 
	- Specifying streams as diagnostic (`aas_garbagecollection` will not touch them)
    - aaworker
      - More efficient folder handling: 
        - aaworker folders are created in ~/.aa (hidden) folder and 
        - aaworker folders are cleaned up regularly after specified time (`aap.options.aaworkercleanup`)
      - aaworker structure is passed during `qsub`

  - Summaries
    - Add default message for "Motion correction summary"
    - Rename "First level contrasts" to "First level thresholded maps"

#### Convenience/Ease of access ####
  - Code development
    - Example module template with session domain ([aamod_session](https://github.com/rhodricusack/automaticanalysis/tree/v4-devel/examples/aamod_session.m))

  - Setting up
    - `aaClass` provides links to ["the aa manuscript"](http://dx.doi.org/10.3389/fninf.2014.00090), to ["the aa website"](http://automaticanalysis.org) and to examples
    - Allow specification of subject in the whole UMS in the same way: `aas_addcontrast`, `aas_addcovariate`, `aas_addevent`, `aamod_dartel_createtemplate`
    - `aap.options.checktasksettingconsistency`: Check whether settings have changed since the last execution and re-run the task accordingly regardless of the doneflag (experimental!)
    - <aa version> and <aa path> are saved in the _aap\_parameters.mat_ file as fields `aap.internal.aapversion` and `aap.internal.aappath`

  - Running
    - `qsub` error provides links to code
    - GUI for `qsub`
    - `aap.options.verbosity`:
      <br>0 - Does not show any message (still crashes in case of error)
      <br>1 - Shows errors only
      <br>2 - Shows every message (default)

  - Results
    - Vertical workflow chart
    - More consistent (and nicer :sparkles:) graph of contrasts 
    - `aamod_firstlevel_threshold` creates overlays along all three axes

## 4.2 ([branch devel-share initial](https://github.com/rhodricusack/automaticanalysis/tree/v4-stable)) ##

Described in ["the aa manuscript"](http://dx.doi.org/10.3389/fninf.2014.00090)

