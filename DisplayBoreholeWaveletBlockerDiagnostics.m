function funHandle = DisplayBoreholeWaveletBlockerDiagnostics(inputData, figHandle)
%DisplayBoreholeWaveletBlockerDiagnostics. Displays the outputs of the blocking algorithm.
%   funHandle = DisplayBoreholeWaveletBlockerDiagnostics(inputData).
%   Displays a figure with 2 subplots.  The top subplot shows the layer
%   importances as determined by taking the average value of the transform
%   in that layer.  This is based on the regions determined by the wavelet
%   blocker algorithm, with regions created by the label function.  The
%   regional averages are calculated by taking the average of every pixel
%   in the region (absolute, so that negative deflections are as important
%   as positive ones).  Layer averages are normalised to 1, with 1 being
%   the most important layer.  Subplot 2 is the same layers plotted with
%   their layer thickness.  This can be a misleading graph, since some
%   layers contain smaller layers (as seen in the wavelet blocking figure).
%
%   This figure is useful when determining which criteria to use to block
%   the bore trace.
%
%
%   Company: Commonwealth Scientific and Industrial Research Organisation
%   (CSIRO), Earth Science and Resource Engineering, 2013
%   Author: Aaron C Davis
%
%   This software is licenced under the Creative Commons Attribution
%   (CC-BY) 3.0 licence (http://creativecommons.org/licenses/by/3.0/)

%% Size the figure
PPI = get(0, 'screenPixelsPerInch');
pWidth = 3.3125;
pHeight = 5;
height = pHeight*PPI;
width = pWidth*PPI;

%% Plot the data
if nargin == 2
  funHandle = figHandle;
  figure(funHandle);
else
  funHandle = figure;
end

clf;
set(funHandle, 'Position', [1*PPI 2*PPI width height]);
set(funHandle, 'PaperUnits', 'inches');
set(funHandle, 'PaperSize', [pWidth pHeight]);
set(funHandle, 'PaperPosition', [0 0 pWidth pHeight]);

aHandle(1) = subplot(2,1,1);
aHandle(2) = subplot(2,1,2);

set(gcf, 'CurrentAxes', aHandle(1));
loglog(inputData.layerMean, 'k.-');
hold on;
% yLim = ylim;
% plot(ones(2,1).*(inputData.nLayer./100), yLim, 'r');

ylabel('Layer importance');
xlabel('Layer number');
title('Layer importance');

set(gcf, 'CurrentAxes', aHandle(2));
loglog(inputData.layerThickness, 'k.-');
ylabel('Layer thickness (m)');
xlabel('Layer number');
title('Layer thickness');

set(funHandle, 'Renderer', 'Painters');
set(funHandle, 'Color', 'w');
