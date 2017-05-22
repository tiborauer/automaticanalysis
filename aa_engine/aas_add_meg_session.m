function [aap]=aas_add_meg_session(aap,name)

% Blank template for a session entry
thissess.name=name;

% And put into acq_details, replacing a single blank entry if it exists
if (length(aap.acq_details.meg_sessions)==1 && isempty(aap.acq_details.meg_sessions.name))
    aap.acq_details.meg_sessions=thissess;
else
    doAdd = true;
    for iSess = 1:numel(aap.acq_details.meg_sessions)
        if strcmp(aap.acq_details.meg_sessions(iSess).name,thissess.name)
            doAdd = false;
        end
    end
    if doAdd, aap.acq_details.meg_sessions(end+1) = thissess; end    
end