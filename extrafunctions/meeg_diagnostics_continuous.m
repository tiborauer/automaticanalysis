function meeg_diagnostics_continuous(EEG,diag,figtitle,savepath)
if nargin < 3, figtitle = 'Sample'; end
if ~isempty(diag.freqrange)
    f = figure('Name',figtitle);
    if isempty(diag.freq)
        pop_spectopo(EEG,1,[],'EEG','freqrange',diag.freqrange,'electrodes','off');
    else
        pop_spectopo(EEG,1,[],'EEG','freq',diag.freq,'freqrange',diag.freqrange,'electrodes','off');
    end
    % title to the main axes
    all_axes = findall(gcf, 'Type', 'axes');
    ax = all_axes(arrayfun(@(a) isequal(a.XLim, diag.freqrange), all_axes));
    title(ax,figtitle);
    if nargin >= 4
        set(f,'Position',[1 1 1920 1080]);
        set(f,'PaperPositionMode','auto');
        print(f,'-noui',savepath,'-djpeg','-r150');
        close(f);
    end
end