function [s w]=aas_runfslcommand(aap,fslcmd)

% Setting paths now done in aa_doprocessing
fslsetup=deblank(aap.directory_conventions.fslsetup);

if not(isempty(fslsetup)) && not(fslsetup(end)==';')
    fslsetup=[fslsetup ';'];
end;

ENV = {...
    'MATLABPATH', path;...
    'FSLOUTPUTTYPE', aap.directory_conventions.fsloutputtype ...
    };

switch (aap.directory_conventions.fslshell)
    case 'none'
        for e = 1:size(ENV,1)
            setenv(ENV{e,1},ENV{e,2});
        end
        cmd=[fslsetup fslcmd];
        [s w]=aas_shell(cmd);
    case {'csh' 'tcsh'}
        cmd=[aap.directory_conventions.fslshell ' -c "'];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('setenv %s %s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd fslcmd '"'];
        [s w]=aas_shell(cmd);
    case 'bash'
        cmd=[aap.directory_conventions.fslshell ' -c "'];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd fslcmd '"'];
        [s w]=aas_shell(cmd);
end;

% Display error if there was one
if (s)
    aas_log(aap,false,sprintf('Error running %s, which was %s',cmd,w));
end;