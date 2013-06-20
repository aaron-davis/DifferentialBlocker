function funHandle = DisplayBoreholeWaveletBlockerLayerSelection(blockedData, dataSwitch, figHandle)
%DisplayBoreholeWaveletBlockerTransform. Display the blocked trace.
%   funHandle = DisplayBoreholeWaveletBlockerLayerSelection(blockedData).
%   Displays the original trace of the data and the result of the wavelet
%   blocker based on the output of the function
%   BoreholeWaveletBlockerLayerSelection.m.  This is the final result of
%   the blocked trace, based on user-specified blocking criteria.  The
%   output is the mean values of the layers that are blocked
%
%   funHandle = DisplayBoreholeWaveletBlockerLayerSelection(blockedData, dataSwitch).
%   Same as above with the additional selection of dataSwitch, a flag
%   operator. dataSwitch causes 2 separate outputs:
%     dataSwitch = 0 (default setting when no other input)
%       Layers are blocked according to the mean of the layer.
%     dataSwitch = 1
%       Layers are blocked according to the median of the layer.
%
%   Company: Commonwealth Scientific and Industrial Research Organisation
%   (CSIRO), Earth Science and Resource Engineering, 2013
%   Author: Aaron C Davis
%
%   This software is licenced under the Creative Commons Attribution
%   (CC-BY) 3.0 licence (http://creativecommons.org/licenses/by/3.0/)

%% Detect the dataSwitch
if nargin < 2
  dataSwitch = 0;
end
switch dataSwitch
  case 0
    fieldName{1} = 'plotLayerMean';
  otherwise
    fieldName{1} = 'plotLayerMedian';
end % dataSwitch

%% Size the figure
PPI = get(0, 'screenPixelsPerInch');
pWidth = 5;
pHeight = 7;
height = pHeight*PPI;
width = pWidth*PPI;

%% Plot the figure
if nargin == 3
  funHandle = figHandle;
  figure(funHandle);
else
  funHandle = figure;
end
clf;
set(funHandle, 'Position', [0.5*PPI 0.5*PPI width height]);
set(funHandle, 'PaperUnits', 'inches');
set(funHandle, 'PaperSize', [pWidth pHeight]);
set(funHandle, 'PaperPosition', [0 0 pWidth pHeight]);

plot(blockedData.data, blockedData.depth, 'k');
hold on;
plot(blockedData.(fieldName{1}), blockedData.plotLayerDepth, 'r');

% Annotation
ylabel('Depth (m)');
xlabel('Data value');
title(['Blocked log (Number of layers = ' num2str(blockedData.nLayer) ')']);
set(funHandle, 'Renderer', 'Painters');
set(funHandle, 'Color', 'w');
set(gca, 'yDir', 'reverse');

grid off;
