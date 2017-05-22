% AA initialisation module - Identify dicom headers from Tim Trio
% Rhodri Cusack MRC CBU Cambridge Dec 2005
% Jul 2012 - can be modified to examine GE or Siemens headers

function  [aap,resp]=aamod_autoidentifyseries(aap,task,i)

resp='';

switch task
    case 'summary'
        resp=[];
        numresp=0;
        if (aap.options.autoidentifyfieldmaps)
            resp{1}='fieldmaps';
            numresp=numresp+1;
        end
        if (aap.options.autoidentifystructural)
            resp{numresp}='structural';
            numresp=numresp+1;
        end
        if (aap.options.autoidentifytmaps)
            resp{numresp}='realtime t-maps';
            numresp=numresp+1;
        end
        switch (numresp)
            case 0
                resp=['No automatic identification done.'];
            case 1
                resp=['Automatically identified ' resp{1} '\n'];
            case 2
                resp=['Automatically identified ' resp{1} ' and ' resp{2} '\n'];
            case 3
                resp=['Automatically identified ' resp{1} ', ' resp{2} ' and ' resp{3} '\n'];
        end
    case 'report'
        
    case 'doit'
        global aaworker
        dontrescan=false;
        aisfn=fullfile(aas_getsubjpath(aap,i),'autoidentifyseries_saved.mat');
        
        % Get a listing of all of the series for this subject
        rawdata_subj=fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).mriname);
        switch (aap.directory_conventions.remotefilesystem)
            case 'none'
                % Looks for DICOMs in here and in all subdirectories
                allpths=genpath(rawdata_subj);
                serieslist=[];
                order_by=[];
                slicelocation=[];
                alldicomfiles={};
                while (length(allpths)>0)
                    [thispth allpths]=strtok(allpths,':');
                    if (isempty(thispth))
                        break;
                    end
                    fn=dir(fullfile(thispth,aap.directory_conventions.dicomfilter));
                    for fnind=1:length(fn)
                        if (~fn(fnind).isdir)
                            fullfn=fullfile(thispth,fn(fnind).name);
                            H=spm_dicom_headers(fullfn);
                            if (isfield(H{1},'SeriesNumber') && isfield(H{1},'AcquisitionNumber'))
                                serieslist=[serieslist H{1}.SeriesNumber];
                                switch(aap.tasklist.currenttask.settings.orderdicomsby)
                                    case 'AcquisitionNumber'
                                        ob=H{1}.AcquisitionNumber;
                                        
                                    case 'InstanceNumber'
                                        ob=str2double(fn(fnind).name(2:end));
                                        if isfield(H{1},'InstanceNumber')
                                            ob=double(H{1}.InstanceNumber);
                                        else
                                            ob=0;
                                        end;
                                end;
                                if isfield(H{1},'SliceLocation')
                                    sl=double(H{1}.SliceLocation);
                                else
                                    sl=0;
                                end;
                                
                                if (H{1}.SeriesNumber>length(alldicomfiles))
                                    alldicomfiles{H{1}.SeriesNumber}=[];
                                    order_by{H{1}.SeriesNumber}=[];
                                    slicelocation{H{1}.SeriesNumber}=[];
                                end
                                alldicomfiles{H{1}.SeriesNumber}{end+1}=fullfn;
                                order_by{H{1}.SeriesNumber}(end+1)=ob;
                                slicelocation{H{1}.SeriesNumber}(end+1)=sl;
                            end
                        end
                    end
                    
                end
                
                rawdata_allseries=unique(serieslist)
                
                for j=1:length(rawdata_allseries);
                    [junk ind]=sort(order_by{rawdata_allseries(j)});
                    alldicomfiles{rawdata_allseries(j)}=alldicomfiles{rawdata_allseries(j)}(ind);
                    
                    
                    aas_log(aap,false,sprintf('Series %d with %d dicom files',rawdata_allseries(j),length(alldicomfiles{rawdata_allseries(j)})));
                end
            case 's3'
                % Use delimiter to get series names as CommonPrefixes
                [aap s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,rawdata_subj,[],'/');
                rawdata_allseries={s3resp.CommonPrefixes.Prefix};
        end
        
        series_spgr=[];
        series_newfieldmap=[];
        series_tmaps=[];
        
        filenumber=0;
        
        tmpdir=aas_gettempfilename();
        
        protocolnames=[];
        
        % Go through each series, and examine type
        for j=1:length(rawdata_allseries)
            % Get the path to a single dicom file from series "j", downloading from S3 first if necessary
            switch(aap.directory_conventions.remotefilesystem)
                case 's3'
                    seriespth=rawdata_allseries{j};
                    while (seriespth(end)==filesep)
                        seriespth=seriespth(1:(end-1));
                    end
                    [pth nme ext]=fileparts(seriespth);
                    dicomseriesname=[nme ext];
                    searchpath=fullfile(rawdata_subj,dicomseriesname);
                    aas_log(aap,false,sprintf('Checking here on s3 %s',searchpath));
                    [aap s3resp]=s3_list_objects(aap,aaworker.bucketfordicom,searchpath,[],[],1);
                    
                    [keypth nme ext]=fileparts(s3resp.Contents(1).Key);
                    dicomfilename=[nme ext];
                    s3_copyfrom_filelist(aap,tmpdir,dicomfilename,aaworker.bucketfordicom,keypth);
                    dicomfilepath=fullfile(tmpdir,dicomfilename);
                    % Option to specify series numbers in terms of numbering of
                    % scanner, or ordering of files
                    
                    if (aap.directory_conventions.rawseries_usefileorder)
                        seriesnum=j;
                    else
                        seriesnum=aas_getseriesnumber(aap,dicomseriesname);
                    end
                case 'none'
                    seriesnum=rawdata_allseries(j);
                    dicomfilepath=alldicomfiles{seriesnum}{1};
            end
            
            % ignoreseries
            if any(aap.acq_details.subjects(i).ignoreseries{d} == seriesnum), continue; end
            
            % For this series, find type from a single DICOM file
            if (~isempty(dicomfilepath))
                hdr=spm_dicom_headers(dicomfilepath);
                
                % Decide whether to ignore this series [djm 20/3/06]
                if any(aap.acq_details.subjects(i).ignoreseries==j)
                    aas_log(aap,false,sprintf('Ignoring series %d',j));
                else
                    aas_log(aap,false,sprintf('Protocol field is %s and got %s',aap.tasklist.currenttask.settings.dicom_protocol_field,hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field)));

                    if (aap.options.autoidentifyfieldmaps)
                        if (findstr(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),aap.directory_conventions.protocol_fieldmap))
                            series_newfieldmap=[series_newfieldmap seriesnum];
                        end
                    end
                    
                    if (aap.options.autoidentifystructural)
                        if (findstr(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),aap.directory_conventions.protocol_structural))
                            if (series_spgr & ...
                                    aap.options.autoidentifystructural_chooselast==0  & ...
                                    aap.options.autoidentifystructural_average==0 & ...
                                    aap.options.autoidentifystructural_multiple==0) %[AVG] for MP2RAGE, for instance
                                aas_log(aap,1,'Automatic series id failed - more than one MPRAGE acquisition was found.');
                            end
                            series_spgr=[series_spgr seriesnum];
                        end
                    end
                    
                    
                    if (aap.options.autoidentifytmaps)
                        % Use directory name rather than protocol to
                        % recognise t maps
                        if (findstr(hdr{1}.(aap.tasklist.currenttask.settings.dicom_protocol_field),'EvaSeries_tTest'))
                            series_tmaps=[series_tmaps seriesnum];
                        end
                    end
                end
            end
            
        end
        % Save file
        [ais_p ais_f ais_e]=fileparts(aisfn);
        aas_makedir(aap,ais_p);
        aapoptions=aap.options;
        if (exist('alldicomfiles','var'))
            save(aisfn,'series_newfieldmap','series_spgr','series_tmaps','aapoptions','alldicomfiles','rawdata_allseries');
        else
            save(aisfn,'series_newfieldmap','series_spgr','series_tmaps','aapoptions');
        end
        
        % Make comment
        comment=[];
        if (aap.options.autoidentifyfieldmaps)
            aap.acq_details.subjects(i).siemensfieldmap={};
            if (length(series_newfieldmap)>=aap.options.autoidentifyfieldmaps_number)
                comment=[comment ' ' sprintf('gre_fieldmapping found %d',series_newfieldmap(1))];
                % Generalisation of fieldmap number...
                for n = 2:aap.options.autoidentifyfieldmaps_number
                    comment=[comment ' ' sprintf(' and %d',series_newfieldmap(n))];
                end
                aap.acq_details.subjects(i).siemensfieldmap=series_newfieldmap;
            else
                aas_log(aap,1,'Automatic series id failed - one of the fieldmap acquisitions was not found.');
            end
        end
        
        
        if (aap.options.autoidentifystructural)
            if (aap.options.autoidentifystructural_chooselast & length(series_spgr)>1)
                series_spgr=series_spgr(length(series_spgr));
            end
            aap.acq_details.subjects(i).structural=series_spgr;
            comment=[comment sprintf(' Structural series %d ',series_spgr)];
        end
        
        if (aap.options.autoidentifytmaps)
            aap.acq_details.subjects(i).tmaps=series_tmaps;
            comment=[comment [' T maps series ' sprintf('%d\t',series_tmaps)]];
        end
        if (length(comment)>0) aas_log(aap,0,comment); end
        
        aap=aas_desc_outputs(aap,i,'autoidentifyseries',aisfn);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
