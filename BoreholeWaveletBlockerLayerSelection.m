function layerData = BoreholeWaveletBlockerLayerSelection(inputData, ...
  layerSelection, layerSelectFlag)
%BoreholeWaveletBlockerLayerSelection. Performs the layer blocking.
%   BoreholeWaveletBlockerLayerSelection(inputData, layerSelection)
%   Blocks the trace based on the layerSelection criteria.  When nargin is
%   2, the layerSelection criteria is a wavelet width value.  The wavelet
%   width is defined as the positive deflection of the wavelet used to
%   determine layer boundaries.  When the user inputs layerSelection (say,
%   2 m), the program searches the wavelet transform for layers that are
%   detectable with a maximum (positive) wavelet width of 2 m.  Then, all
%   layers that are under that maximum size get traced.  The layer boundary
%   is found by tracing the layer regional contour back to the LHS of the
%   wavelet matrix.  Layers are selected based on this.
%
%   BoreholeWaveletBlockerLayerSelection(inputData, layerSelection, layerSelectFlag)
%   Same as above, with the addition of a layerSelectFlag.  Currently,
%   there are 5 options:
%     layerSelectFlag = 0 (default)
%       Layers are selected on a maximum wavelet width as input in
%       layerSelection.
%     layerSelectFlag = 1 
%       The program looks for first 1:layerSelection layers, based on
%       importance, where layerSelection should be an integer.  For
%       example, if layerSelectFlag == 1, and layerSelection = 10, the program
%       will select the 10 most important layers in the trace, which are
%       definied by the average deflection of the layer from the wavelet
%       transform.
%     layerSelectFlag = 2
%       The program will select the first layerSelection layers, where
%       layerSelection is taken as a percentage of all layers found by the
%       blocking algorithm.  For example, if layerSelectFlag == 2, and 
%       layerSelection = 10, the program will look for the top 10% of
%       important layers.  If the data is very noisy, there will be a lot
%       of layers, so care must be taken when using this selection.
%     layerSelectFlag = 3
%       The program will select all layers that are greater than or equal
%       to the specified layer thickness given by the value layerSelection.
%     layerSelectFlag = 4
%       The program will select all layers that are less than or equal
%       to the maximum specified layer thickness given by the value layerSelection.
%
%   NOTE:
%   Despite the layer selection criteria looking for layers based on
%   wavelet width, # of important layers, etc, there may be some layers of
%   the trace that are not covered by this selection.  In such an example,
%   the remiander of the log trace is blocked off, creating extra layers so
%   that we can fill the entire trace.  Let's propose that we have a trace
%   that is 20 m deep, and there is a large delfection at 10-12m, creating
%   a very important layer.  If we selected layerSelectFlag == 1 and
%   layerSelection = 1, meaning that we want to locate only the most
%   important layer, we will need to artificially block of 0-10m and
%   12-20m, thereby creating 3 layers.  
%
%   layerData = struct(...
%     'depth', inputData.depth, ... Depth of the trace, repeated
%     'data', inputData.data, ... Data from the trace, repeated
%     'nLayer', nLayers, ... # of layers found by blocking
%     'layerSelection', layerSelectFlag, ... layer selection (width, # or %)
%     'layerDepth', layerDepth, ... depths of the layers
%     'layerThickness', layerThickness, ... thickness of layers
%     'layerMean', layerMean, ... mean vlaue of the trace for each layer
%     'layerMedian', layerMedian, ... median vlaue of the trace for each layer
%     'layerVariance, layerVariance ... variance of the trace for each layer
%     'plotLayerDepth', plotDepth, ... plot-friendly layer depth
%     'plotLayerMean', plotDataMean, ... plot-friendly layer mean
%     'plotLayerMedian', plotDataMedian) plot-friendly layer median
%
%   Company: Commonwealth Scientific and Industrial Research Organisation
%   (CSIRO), Earth Science and Resource Engineering, 2013
%   Author: Aaron C Davis
%
%   This software is licenced under the Creative Commons Attribution
%   (CC-BY) 3.0 licence (http://creativecommons.org/licenses/by/3.0/)

%% Determine the layerSelectFlag
if nargin < 3
  layerSelectFlag = 0;
end

%% Select layers based on layerSelectFlag
switch layerSelectFlag
  case 0 % We select layers based on a wavelet width
    waveletNumber = min(find(inputData.waveletWidth<=layerSelection, 1, 'last'), ...
      inputData.nWavelet);
    layersInIndex = unique(inputData.layerLabel(:, waveletNumber));
  case 1 % We select the first layerSection number of important layers
    layersInIndex = 1:floor(layerSelection);
  case 2% We select the first layerSelection percentage of layers
    layersInIndex = 1:floor(layerSelection.*inputData.nLayer./100);
  case 3
    layersInIndex = find(inputData.layerThickness >= layerSelection);
  otherwise 
    layersInIndex = find(inputData.layerThickness <= layerSelection);
    
end % switch layerSelectionFlag

% These are the indexes for the layers that we are going to use in the
% blocker
layerDepthIndexSorted = inputData.layerDepthIndex(layersInIndex, :);

%% Find the indexes and the depths of the blocking
% We have to sort the layers and make them continuous.  Some layers are
% contained in other more important layers.  We need to create extra
% layering to allow for this to occur, and to make it continuous.
allLayers = layerDepthIndexSorted(:);
blockedLayers = [];
[findNextLayerBoundary nextIndex] = min(allLayers);
blockedLayers(1, 1) = findNextLayerBoundary;
allLayers(nextIndex) = [];
[findNextLayerBoundary nextIndex] = min(allLayers);
blockedLayers(1, 2) = findNextLayerBoundary;
allLayers(nextIndex) = [];
blockedLayers(2, 1) = blockedLayers(1, 2) +1;
layerIndex = 2;
while ~isempty(allLayers)
  [findNextLayerBoundary nextIndex] = min(allLayers);
  if findNextLayerBoundary == blockedLayers(layerIndex, 1)
      allLayers(nextIndex) = [];  
      [findNextLayerBoundary nextIndex] = min(allLayers);
      blockedLayers(layerIndex, 2) = findNextLayerBoundary;
  else
    blockedLayers(layerIndex, 2) = findNextLayerBoundary;
  end % if
  allLayers(nextIndex) = [];  
  layerIndex = layerIndex + 1;
  blockedLayers(layerIndex, 1) = blockedLayers(layerIndex-1, 2) + 1;  
end % iLayer
blockedLayers(end,:) = [];
% Check if the blocking algorithm starts from the beginning index, and end
% at the final index to ensure the entire trace is covered by the blocking
% program,
if blockedLayers(1,1) ~= 1
  blockedLayers = [1 blockedLayers(1,1)-1; blockedLayers];
end %if blockedLayers
if blockedLayers(end,2) ~= inputData.nData
  blockedLayers = [blockedLayers; blockedLayers(end, 2) + 1 inputData.nData];
end %if blockedLayers

% Calculate the depths and thickness of the layers found
nLayers = size(blockedLayers, 1);
layerDepth = inputData.depth([blockedLayers(:,1); blockedLayers(end,2)]);
layerThickness = diff(layerDepth);

%% Blocking data
% We block the data according to the mean and median of the layer values.
% Also output is the data in a plotting-friendly manner, so that the
% blocked layers appear as a stair-plot.
plotDataMean = [];
plotDataMedian = [];
plotDepth = [];
layerMean = zeros(nLayers, 1);
layerMedian = zeros(nLayers, 1);

for iLayer = 1:size(blockedLayers, 1)
  layerIndex = blockedLayers(iLayer,1):blockedLayers(iLayer,2);
  layerMean(iLayer, 1) = mean(inputData.data(layerIndex));
  layerMedian(iLayer, 1) = median(inputData.data(layerIndex));
  layerVar(iLayer, 1) = var(inputData.data(layerIndex));
  plotDataMean = [plotDataMean; repmat(layerMean(iLayer, 1), 2, 1)];
  plotDataMedian = [plotDataMedian; repmat(layerMedian(iLayer, 1), 2, 1)];
  plotDepth = [plotDepth; repmat(layerDepth(iLayer), 2, 1)];
end % iLayer

plotDepth(1) = [];
plotDepth(end+1) = layerDepth(end);

%% Output the data in a struct
layerData = struct(...
  'depth', inputData.depth, ...
  'data', inputData.data, ...
  'nLayer', nLayers, ...
  'layerSelection', layerSelectFlag, ...
  'layerDepth', layerDepth, ...
  'layerThickness', layerThickness, ...
  'layerMean', layerMean, ...
  'layerMedian', layerMedian, ...
  'layerVariance', layerVar, ...
  'plotLayerDepth', plotDepth, ...
  'plotLayerMean', plotDataMean, ...
  'plotLayerMedian', plotDataMedian);
  