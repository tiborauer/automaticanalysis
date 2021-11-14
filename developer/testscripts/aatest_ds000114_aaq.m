function aatest_ds000114_aaq(parameterfile, deleteprevious, wheretoprocess)

% developer PR test script
%
% description: test aaq (batch)
% dataset: ds000114

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

% -------------------------------------------------------------------------
% analysis options - batch
% -------------------------------------------------------------------------

aap.options.wheretoprocess = 'batch'; % localsingle
aap.directory_conventions.poolprofile = 'local';
aap.options.aaparallel.numberofworkers = 2;
aap.options.aaparallel.memory = 2;
aap.options.aaparallel.walltime = 1;
aap.options.aaworkermaximumretry = 0; % do not retry; it would remove files
aap.options.aaworkerGUI = 0;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

% anatomy data is in session 3
aap = aas_processBIDS(aap,[],{'finger_foot_lips','line_bisection'},{'sub-01','sub-02'});

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);

