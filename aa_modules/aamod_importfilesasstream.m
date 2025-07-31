% Import files as stream
% This module is typically added automatically by aas_addinitialstream  in
% your user script

function [aap,resp]=aamod_importfilesasstream(aap,task,varargin)

resp='';

switch task
    case 'report'
    case 'doit'
        streamname=aap.tasklist.currenttask.outputstreams.stream{1};
        
        foundmatch=false;
        matches=aap.tasklist.currenttask.settings.match;
        for matchind=1:length(matches)
            % Find where this module is working
            fti=matches(matchind).filenames;
            switch(length(varargin))
                case 0
                    foundmatch=true;
                    currpth=aas_getstudypath(aap);
                case 1
                    if (varargin{1}==matches(matchind).subject)
                        foundmatch=true;
                        currpth=aas_getsubjpath(aap,varargin{1});
                        break;
                    end
                case 2
                    if (varargin{1}==matches(matchind).subject && varargin{2}==matches(matchind).session)
                        foundmatch=true;
                        currpth=aas_getsesspath(aap,varargin{1},varargin{2});
                        break;
                    end
            end
        end
        
        if (~foundmatch)
            aas_log(aap,false,'Warning: No match for import - do you have right number of aas_addinitialstream commands?');
            return;
        else
            aas_log(aap,false,sprintf('Importing %d files',length(fti)));
        end
        % Copy the files to import into the current module's directory
        allnames={};
        localnames={};
        for fileind=1:length(fti)
            ext = spm_file(fti{fileind}, 'ext');
            nme = spm_file(fti{fileind}, 'basename');
            if ~isempty(ext), ext = ['.' ext]; end
            while ~isempty(spm_file(nme, 'ext'))
                newext = spm_file(nme, 'ext');
                nme = spm_file(nme, 'basename');
                ext = ['.' newext ext];
            end
            timesalready=sum([strcmp([nme ext],allnames)]);
            allnames{end+1}=[nme ext];
            localnames{fileind}=fullfile(currpth,sprintf('%s-%d%s',nme,timesalready,ext));
            if aas_getsetting(aap,'unzip')
                switch aas_getsetting(aap,'unzip')
                    case 'gunzip'
                        localnames{fileind} = strrep(localnames{fileind},'.gz','');
                        aas_shell(sprintf('gunzip -dc %s > %s',fti{fileind},localnames{fileind}));
                end
            else
                copyfile(fti{fileind},localnames{fileind});
            end
        end
        
        % Describe the output stream
        switch(length(varargin))
            case 0
                aap=aas_desc_outputs(aap,streamname,localnames);
            case 1
                aap=aas_desc_outputs(aap,varargin{1},streamname,localnames);
            case 2
                aap=aas_desc_outputs(aap,varargin{1},varargin{2},streamname,localnames);
        end
        
    case 'checkrequirements'
        streamname=aap.tasklist.currenttask.outputstreams.stream{1};
        matches=aap.tasklist.currenttask.settings.match;
        for matchind=1:length(matches)
            
            fti=matches(matchind).filenames;
            if (iscell(fti))
                for fileind=1:length(fti)
                    if (~exist(fti{fileind},'file'))
                        aas_log(aap,true,sprintf('File marked as initial import for stream %s not found: %s',streamname,fti{fileind}));
                    end
                end
            else
                for fileind=1:size(fti,1)
                    if (~exist(fti(fileind,:),'file'))
                        aas_log(aap,true,sprintf('File marked as initial import for stream %s not found: %s',streamname,fti(fileind,:)));
                    end
                end
            end
        end
    otherwise
        aas_log(aap, 0, sprintf('%s: Unknown task %s',mfilename, task));
end














