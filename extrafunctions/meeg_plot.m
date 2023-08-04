function h = meeg_plot(cfg,data)
% cfg
%   - parameter - char
%   - channels  - {1xN}
%   - latency   - [Nx2], second
%   - view      - char (item from VAs
%   - statlim   - [Nx2]
%   - scaling   - 'joint','individual'
%   - showlabel - true/false

%% Contants
FIGWIDTH = 1080;
if isfield(cfg,'width'), FIGWIDTH = cfg.width; end
TOPO_MAXCOL = 4; % maximum number of columns of plot mosaic, should be even

VAs = {...
    'A'   { 180, 0};...
    'P'   {   0, 0};...
    'S'   {   0, 90};...
    'I'   { 180,-90};...
    'R'   {  90,  0};...
    'L'   { -90,  0};...
    'RAS' { 135, 45};...
    'LAS' {-135, 45};...
    'RPS' {  45, 45};...
    'LPS' { -45, 45};...
    'RAI' { 135,-45};...
    'LAI' {-135,-45};...
    'RPI' {  45,-45};...
    'LPI' { -45,-45};...
    };

%% Groups
cfg.showlabel = ~isfield(cfg,'showlabel') || cfg.showlabel;
if numel(data) > 1
    labels = arrayfun(@(x) sprintf('group #%d',x),1:numel(data),'UniformOutput',false);
else
    labels = {'group'};
end


%% Stats
if isfield(data{1},'stat')
    stat = data{1}.stat;
    data{1} = rmfield(data{1},'stat');
    
    stat = keepfields(stat,{'stat' 'mask' 'cfg' 'dimord' 'posclusters' 'posclusterslabelmat' 'negclusters' 'negclusterslabelmat'});
    if ndims(stat.stat) ~= ndims(data{1}.(cfg.parameter)) % assume only cases with singleton dimension in stat        
        statdims = size(stat.stat);
        statdims(statdims==1) = [];
        if numel(statdims) ~= ndims(data{1}.(cfg.parameter)) || ~all(statdims == size(data{1}.(cfg.parameter)))            
            aas_log([],true,['extra non-singleton dimension "' statdims{d} '" found in stat'])
        end
        for f = fieldnames(stat)'
            stat.(f{1}) = squeeze(stat.(f{1}));
        end
    end
    
%     % TODO - make use of separate positive and negative clusters
%     if strcmp(stat.cfg.correctm,'cluster')
%         if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
%             pos_cluster_pvals = [stat.posclusters(:).prob];
%             pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
%             pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
%         else
%             pos = stat.mask*0;
%         end
%         
%         if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
%             neg_cluster_pvals = [stat.negclusters(:).prob];
%             neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
%             neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
%         else
%             neg = stat.mask*0;
%         end
%     else
%         pos = stat.mask & (stat.stat > 0);
%         neg = stat.mask & (stat.stat < 0);
%     end
   
    data = [data(1) data];
    data{1}.(cfg.parameter) = stat.stat;
    data{1}.mask = logical(stat.mask);
    for s = 2:numel(data)
        if isfield(data{s},'mask'), data{s}.mask = logical(stat.mask); end
    end
    labels = [{'stat'} labels];
end
TOPO_MAXCOL = floor(TOPO_MAXCOL/numel(data));

%% Time
if isfield(data{1},'time')
    if isempty(cfg.latency), cfg.latency = [data{1}.time(1) data{1}.time(end)]; end
    nlPlot = ceil(size(cfg.latency,1)/TOPO_MAXCOL);
    sizeplot = [nlPlot size(cfg.latency,1)*numel(data)/nlPlot];
else
    cfg.latency = [0 0];
    for s = 1:numel(data)
        data{s}.time = 0;
    end
    sizeplot = [1,numel(data)];
end

mplier = 1;
if isfield(cfg,'channels'), mplier = numel(cfg.channels); end
if isfield(cfg,'view') && strcmp(cfg.view,'ortho'), mplier = 3; end
sizeplot(1) = sizeplot(1)*mplier;

%% Scaling
for s = 1:numel(data)
    if isfield(data{s},'mask') && any(data{s}.mask(:))
        mask = data{s}.mask;
    else
        mask = true(size(data{s}.(cfg.parameter)));
    end
    if strcmp(labels{s},'stat')
        if isfield(cfg,'statlim') && ~isempty(cfg.statlim)
            minval(s) = cfg.statlim(1);
            maxval(s) = cfg.statlim(2);
        else
            minval(s) = prctile(data{s}.(cfg.parameter)(mask),1,'all');
            maxval(s) = prctile(data{s}.(cfg.parameter)(mask),99,'all');
        end
    else
        if isfield(cfg,'scaling') && strcmp(cfg.scaling,'individual')
            minval(s) = prctile(data{s}.(cfg.parameter)(mask),1,'all');
            maxval(s) = prctile(data{s}.(cfg.parameter)(mask),99,'all');
        else
            minval(s) = prctile(cell2mat(cellfun(@(d) d.(cfg.parameter)(mask), data(~strcmp(labels,'stat')), 'UniformOutput',false)),1,'all');
            maxval(s) = prctile(cell2mat(cellfun(@(d) d.(cfg.parameter)(mask), data(~strcmp(labels,'stat')), 'UniformOutput',false)),99,'all');
        end
    end
    if minval(s) == maxval(s)
        minval(s) = minval(s)-abs(0.1*minval(s)); 
        maxval(s) = maxval(s)+abs(0.1*maxval(s));
    end
    % colormaps
    if (minval(s) < 0) && (maxval(s) > 0)
        r = maxval(s)/-minval(s);
        if r > 1, cmaps{s} = [cmapcold(round(64/r)); cmaphot(64)];
        else, cmaps{s} = [cmapcold(64); cmaphot(round(r*64))];
        end
    elseif minval(s) < 0, cmaps{s} = cmapcold(64);
    else
        cmaps{s} = cmaphot(64);
    end
end

%% Plot
h = figure;
set(h, 'color', [1 1 1]);
for t = 1:size(cfg.latency,1)
    for s = 1:numel(data)
        indAx = (t-1)*numel(data)+s;
        indAx = rem(indAx-1,sizeplot(2))+1 + floor((indAx-1)/sizeplot(2))*mplier*sizeplot(2);
        p = subplot(sizeplot(1),sizeplot(2),indAx, 'Parent',h);
        
        tmpcfg = [];
        tmpcfg.latency = cfg.latency(t,:);
        if isnumeric(cfg.latency)
            strTitle = sprintf('%03d-%03d ms',cfg.latency(t,:)*1000);
            tmpcfg.avgovertime = 'yes';
        else % all
            strTitle = 'all';
        end
        if ~isempty(labels{s}), strTitle = sprintf('%s: %s',labels{s},strTitle); end
        dataPlot = ft_selectdata(tmpcfg,rmfield(data{s},intersect(fieldnames(data{s}),{'inside','mask'})));
        if isfield(data{s},'mask'), dataPlot.mask = data{s}.mask; end
        if isnumeric(cfg.latency)
            dataPlot.dimord = strrep(dataPlot.dimord,'_time',''); 
            dataPlot = rmfield(dataPlot,'time');
        end
        
        if isfield(dataPlot,'mask') % stats
            if ~any(dataPlot.mask(:))
                strTitle = sprintf('%s\nns.',strTitle);
            else
                dataPlot.(cfg.parameter) = dataPlot.(cfg.parameter).*dataPlot.mask;
            end
        end
        
        tmpcfg = keepfields(cfg,{'layout' 'parameter'});
        tmpcfg.interactive = 'no';
        if isfield(dataPlot,'label')
            tmpcfg.zlim = sort([minval(s) maxval(s)]);
            if ft_datatype(dataPlot,'freq') % TFR
                switch numel(strsplit(dataPlot.dimord,'_'))
                    case 3
                        dims = strsplit(dataPlot.dimord,'_');
                        if ~strcmp(dims{3},'time')
                            dataPlot.time = dataPlot.(dims{3});
                            dataPlot.dimord = 'chan_freq_time';
                        end                        
                        
                        if ~isfield(cfg,'channels')  % multiplot
                            adjustaxes = true;
                            FIGWIDTH = 1.5*1080;
                            
                            tmpcfg.showlabels = 'yes';
                            tmpcfg.showoutline = 'yes';
                            ft_multiplotTFR(tmpcfg,dataPlot);
                            currfig = gcf;
                            arrayfun(@(x) copyobj(x,p), get(get(currfig,'CurrentAxes'),'Children'));
                        else
                            adjustaxes = false;
                            for ch = 1:numel(cfg.channels)
                                p(ch) = subplot(sizeplot(1),sizeplot(2),indAx + (ch-1)*sizeplot(2), 'Parent',h);
                                tmpcfg.channel = cfg.channels(ch);
                                ft_singleplotTFR(tmpcfg,dataPlot);
                                currfig(ch) = gcf;
                                arrayfun(@(x) copyobj(x,p(ch)), get(get(currfig(ch),'CurrentAxes'),'Children'));
                                ylabel(p(ch),cfg.channels{ch});
                            end
                        end
                    case 2 %multiplot
                        adjustaxes = true;

                        % add time
                        dataPlot.(cfg.parameter)(:,:,2) = dataPlot.(cfg.parameter);
                        dataPlot.time = [0 1];
                        dataPlot.dimord = [dataPlot.dimord '_time'];
                        
                        tmpcfg.xlim = dataPlot.time([1 end]);
                        tmpcfg.comment = 'yes';
                        ft_multiplotTFR(tmpcfg,dataPlot);
                        currfig = gcf;
                        arrayfun(@(x) copyobj(x,p), get(get(currfig,'CurrentAxes'),'Children'));
                    case 1 %topoplot
                        adjustaxes = true;
                        FIGWIDTH = 0.75*1080;
                        
                        tmpcfg.showoutline = 'yes';
                        tmpcfg.marker = 'labels';
%                         tmpcfg.marker = 'off';
                        tmpcfg.markerfontsize = 4;
                        tmpcfg.comment = 'no';
%                         tmpcfg.highlight = 'on';
%                         tmpcfg.highlightchannel = dataPlot.label(data{1}.mask);
%                         tmpcfg.highlightsymbol = '*';
%                         tmpcfg.highlightcolor = [1 1 1];
                        dataPlot = keepfields(dataPlot,{'label',tmpcfg.parameter,'elec','dimord'});
                        figure; ft_topoplotER(tmpcfg,dataPlot);
                        currfig = gcf;
                        arrayfun(@(x) copyobj(x,p), flipud(get(get(currfig,'CurrentAxes'),'Children')));
                end
            elseif ft_datatype(data,'timelock')  % ER
            end
            
        elseif isfield(dataPlot,'labelcmb') % connectivity
            adjustaxes = false;
            if ndims(dataPlot.(cfg.parameter)) == 2, cmaps{s} = []; end
            tmpcfg = keepfields(cfg,'parameter');
            figure; ft_connectivityplot(tmpcfg,dataPlot);
            currfig = gcf;
            arrayfun(@(x) copyobj(x,p), flipud(get(get(currfig,'CurrentAxes'),'Children')));
        
        elseif isfield(dataPlot,'dim')
            tmpcfg = keepfields(cfg,'latency');
            tmpcfg.funparameter = cfg.parameter;
            tmpcfg.funcolormap = cmaps{s};
            tmpcfg.funcolorlim = [minval(s) maxval(s)];
            if isfield(dataPlot,'mask') && any(dataPlot.mask), tmpcfg.maskparameter = 'mask'; end
            switch cfg.view
                case 'ortho'
                    adjustaxes = true;
                    if ~exist('peakposvox','var')
                        peakposvox = data{1}.transform\[dataPlot.pos(find(dataPlot.((cfg.parameter))==max(dataPlot.(cfg.parameter)),1,'first'),:) 1]';
                    end                    
                    tmpcfg.method = 'slice';
                    tmpcfg.nslices = 1;
                    for a = 1:3
                        p(a) = subplot(sizeplot(1),sizeplot(2),indAx + (a-1)*sizeplot(2), 'Parent',h);
                        tmpcfg.slicerange = [peakposvox(a) peakposvox(a)];
                        tmpcfg.slicedim = a;
                        ft_sourceplot(tmpcfg, dataPlot);
                        colormap(cmaps{s});
                        currfig(a) = gcf;
                        arrayfun(@(x) copyobj(x,p(a)), flipud(get(get(currfig(a),'CurrentAxes'),'Children')));
                        if isfield(tmpcfg,'maskparameter'), p(a).Children(1).AlphaData(p(a).Children(1).CData==0) = 0; end
                        if a == 2, view(p(a),[90 -90]); end
                    end
            end
            
        elseif isfield(dataPlot,'tri')
            adjustaxes = true;
            tmpcfg = keepfields(cfg,'latency');
            tmpcfg.funparameter = cfg.parameter;
            tmpcfg.funcolormap = cmaps{s};
            tmpcfg.funcolorlim = [minval(s) maxval(s)];
            tmpcfg.method = 'surface';
            tmpcfg.camlight = 'no';
            ft_sourceplot(tmpcfg, dataPlot);
            currfig = gcf;
            arrayfun(@(x) copyobj(x,p), flipud(get(get(currfig,'CurrentAxes'),'Children')));
            view(p,VAs{strcmp(VAs(:,1),cfg.view),2}{:});
            camlight(p);
        end
        for i = 1:numel(p)
            cellfun(@(lim) set(p(i),lim,get(get(currfig(i),'CurrentAxes'),lim)), {'XLim','YLim','ZLim','CLim'});
            if adjustaxes
                axis(p(i),'image');
                set(p(i), 'Tag', 'ik', 'Visible', 0);
                set(p(i), 'Tag', 'jk', 'Visible', 0);
                set(p(i), 'Tag', 'ij', 'Visible', 0);
            end
            if cfg.showlabel
                title(p(1), strTitle);
                set(get(p(1),'Title'), 'Visible', 1);
            end
            if ~isempty(cmaps{s})
                colormap(p(i),cmaps{s});
                cb = colorbar(p(i));
                cb.FontSize = 16;
                cb.FontWeight = 'bold';
                if isfield(dataPlot,'mask') && any(dataPlot.mask)
                    pow10 = round(log10(max(abs(str2double(cb.TickLabels)))/max(abs(dataPlot.(cfg.parameter)(dataPlot.mask)))));
                else
                    pow10 = round(log10(max(abs(str2double(cb.TickLabels)))/max(abs(dataPlot.(cfg.parameter)))));
                end
                if abs(pow10) > 1
                    cb.TickLabels = cellstr(num2str(str2double(cb.TickLabels)/10^pow10));
                end
            end
            close(currfig(i));
        end        
    end
end
set(h,'Position',[100,100,FIGWIDTH,round(sizeplot(1)/sizeplot(2)*FIGWIDTH)]);
set(h,'PaperPositionMode','auto');

end

%% Colormap functions
function cmap = cmapbase
[s, FT] = aas_cache_get([],'fieldtrip');
if s
    FT.load;
    FT.addExternal('brewermap');
    cmap = flipud(brewermap(128,'RdBu'));
else
    cmap = flipud([create_grad([0.4 0 0.1],[1 0 0],32);create_grad([1 0 0],[1 1 1],32);create_grad([1 1 1],[0 0 1],32);create_grad([0 0 1],[0 0.2 0.4],32)]);
end
end

function cmap = cmaphot(n)
cmap = cmapbase;
cmap = cmap(65:(end-(64-n)),:);
end

function cmap = cmapcold(n)
cmap = cmapbase;
cmap = cmap((64-n+1):64,:);
end