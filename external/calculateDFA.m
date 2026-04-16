function [DFAExponent, meanDF, windowSizes] = calculateDFA(Signal, windowSizes, windowOverlap)
%     Originally created by Richard Hardstone (2020), rhardstone@gmail.com

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Calculates DFA (on a set of window sizes) for signal
% Steps refer to description of DFA algorithm in Figure 1 of: 
    % Detrended fluctuation analysis: a scale-free view on neuronal oscillations 
    % Hardstone, Richard, Simon-Shlomo Poil, Giuseppina Schiavone, Rick Jansen, Vadim V. Nikulin, Huibert D. Mansvelder, and Klaus Linkenkaer-Hansen 
    % Frontiers in physiology 3 (2012): 450.


% Signal dimensions (numSamples,numChannels)
% windowSize in samples
% windowOverlap is fraction of overlap between windows (0-1)
% Signal dimensions (numSamples,numChannels)
% windowSize in samples
% windowOverlap in samples


lengthSignal = size(Signal,1);
numChannels = size(Signal,2);
meanDF = zeros(numChannels,length(windowSizes));
DFAExponent = zeros(numChannels,1);
windowSizes = windowSizes(:);

for i_channel = 1:size(Signal,2)
    for i_windowSize = 1:length(windowSizes)
        windowOffset = floor(windowSizes(i_windowSize) * (1-windowOverlap));
        allWindowIndex = createWindowIndices(lengthSignal, windowSizes(i_windowSize), windowOffset);
        originalAmplitude = Signal(:,i_channel);
        signalProfile = cumsum(originalAmplitude - mean(originalAmplitude));    %Step A->B
        xSignal = signalProfile(allWindowIndex); %Step B->C
        dSignal = detrend(xSignal'); %Step C
        w_detrendedFluctuations = std(dSignal,1); %std(dSignal,1) means normalized by N instead of N-1
        meanDF(i_channel,i_windowSize) = mean(w_detrendedFluctuations); 
    end
   
    X = [ones(length(windowSizes),1) log10(windowSizes)];
    Y = log10(meanDF(i_channel,:))';
    if length(windowSizes) > 1
        regressOutput = regress(Y,X);
        DFAExponent(i_channel) = regressOutput(2,1); %Step D
    else
        DFAExponent(i_channel) = nan;
    end
end

function allWindowIndex = createWindowIndices(lengthSignal, lengthWindow, windowOffset)
%Gets indices for a set of windows of size (lengthWindow) with a set overlap between them

windowStarts = (1:windowOffset:lengthSignal-lengthWindow+1)-1;
numWindows = length(windowStarts);

oneWindowIndex = 1:lengthWindow;
allWindowIndex = repmat(oneWindowIndex,[numWindows,1]);

allWindowIndex = allWindowIndex + repmat(windowStarts',[1 lengthWindow]);

