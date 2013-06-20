function [outputData, layerMeanMatrix, waveletMatrix] = ...
  BoreholeWaveletBlocker(inputData, inputDepth)
%BoreholeWaveletBlocker Divide a borehole trace into regions.
%   outputData = BoreholeWaveletBlocker(inputData, inputDepth)
%   converts the borehole trace, as defined by the inputDepth and the
%   inputData, into a matrix of regions that define layer boundaries based
%   on a set of wavelets that are created using the CreateWaveletDouble
%   function.  These functions are based on a wavelet that starts as 
%   [-2/9 -4/9 0 6/9 6/9 0 -4/9 -2/9] and steps up in 4.  The wavelet starts at 8. 
%   All wavelets would have positive peak of 1, but are even numbered so that
%   the peak of 1 is never really there.  They are, essentially,
%   even-numbered modifications of the primary wavelet [-1/2 1 -1/2], which
%   is a double derivative operator in the simplest form.
%   The output of the wavelet blocking program is a structure,
%   'outputData':
%   'depth' depth values.
%   'data' data values. 
%   'depthStep' the mean depth step. 
%   'nData' the number of data in the array. 
%   'nWavelet' the number of wavelets used to calculate derivatives.
%   'nLayer' the number of layers found by the routine.
%   'waveletTransform' a matrix of the wavelet transform data 
%   (size nData x nWavelet).
%   'layerLabel' a matrix of size nData x nWavelet that contains the layers 
%   labelled according to layer importance (most important being 1).
%   'layerMean' values of the mean of the layers.  This is the average value of the
%   wavelet transform over the entire region calculated from the layerLabel
%   matrix.  Means are calculated based on area, so a small layer with a
%   large change in geophysical value will have large values (based on the
%   deflection).  These values determine the layer importance, since large
%   changes in the values of the trace have significance in terms of
%   defining a layer.
%   'layerDepthIndex' indeces in depth array that signify the beginning and
%   end of layers, (based on importance).
%   'layerDepth' depth of the layer boundaries.
%   'layerThickness' thickness of the layers, based on importance.
%   'waveletWidth' an array that contains the wavelet widths used in the 
%   transform.  The wavelet widths are based on the positive side of the 
%   wavelet used and the depth step value (so that the wavelet width on the 
%   wavelet that is 8 taps is (8+1)./3 = 3, multiplied by the depth step.
%   outputData = struct(...
%     'depth', inputDepth, ...
%     'data', inputData, ...
%     'depthStep', depthStep, ...
%     'nData', nData, ...
%     'nWavelet', nWavelet, ...
%     'nLayer', nLayer, ...
%     'waveletTransform', diffData, ...
%     'layerLabel', layerLabelMatrix, ...
%     'layerMean', layerSort, ...
%     'layerDepthIndex', layerDepthIndexSorted, ...
%     'layerDepth', layerDepth, ...
%     'layerThickness', layerThickness, ...
%     'waveletWidth', waveletDerivativeWidths);
%
%   [outputData layerMeanMatrix] = BoreholeWaveletBlocker(inputData, inputDepth)
%   Also returns the 'outputData' structure, but also a matrix of size
%   nData x nWavelet.  The layerMeanMatrix is the same as the layerLabel
%   matrix output in outputData struct except that the layer labels are
%   replaced with the layerMean value.  This is useful for a display
%   purpose.
%
%   [outputData layerMeanMatrix waveletMatrix] = ...
%                                 BoreholeWaveletBlocker(inputData, inputDepth)
%   Same outputs as above but with a matrix of size 2*nData+2 x nWavelet.
%   This matrix contains all the wavelets used in the transform.
%
%   NOTE:
%   This code relies on an external function: ilabel.  The function ilabel is written by 
%    Perter Corke and is available from following link:
%    http://www.petercorke.com/Machine_Vision_Toolbox.html
%   Please view his licencing for future works.
%
%   Company: Commonwealth Scientific and Industrial Research Organisation
%   (CSIRO), Earth Science and Resource Engineering, 2013
%   Author: Aaron C Davis
%
%   This software (except the fucntion ilabel.m) is licenced under the Creative 
%   Commons Attribution (CC-BY) 3.0 licence (http://creativecommons.org/licenses/by/3.0/)
%   If you intend to use this software for commercial purposes, you may
%   need to change the function ilabel.m for some other sort of CC-BY code.

%% Prepare the dataset
% Determine depth direction
if mean(diff(inputDepth)) < 0
  inputDepth = flipud(inputDepth);
  inputData = flipud(inputData);
end
depthStep = mean(diff(inputDepth));

% Make the data set even numbered
nPoint = numel(inputData);
if mod(nPoint,2) == 1
  nPoint = nPoint - 1;
  inputData(end) = [];
  inputDepth(end) = [];
end % if

%% Prepare the data for wavelet blocking
tData = [0; -flipud(inputData-mean(inputData)); 0; inputData-mean(inputData)];

%% Prepare the FFT of the borehole data
nPointTData = numel(tData);
tDataFFT = fft(tData);

%% Compute the wavelet transform
fprintf(1, '%s', 'Performing the fast wavelet transforms...');
nData = numel(inputData);
waveletWidth = 8:4:min((3.*nData-1), nPointTData);
nWavelet = numel(waveletWidth);
% With the orthogonal wavelet, the width of the positive portion of the hat
% is (hatWidth)./3 wide.  For example, if the total wavelet functional
% length is 8, the positive portion of the hat is (8 + 1)./3 = 3 points wide.
% This can create a slight problem with layer boundaries when we go to
% small numbers.

diffData = zeros(nData, nWavelet);
waveletMatrix = zeros(nPointTData, nWavelet); 

counter = 1;
for iWave = waveletWidth
  waveTrain = CreateWaveletDouble(iWave, nPointTData)';
  waveletMatrix(:, counter) = waveTrain;
  % FFT of the wavelet
  waveTrainFFT = fft(waveTrain, nPointTData);
  % Compute the convolution & transform back
  tDiffData = ifft(waveTrainFFT.*tDataFFT);
  % Take only one bore length worth of data
  diffData(:, counter) = tDiffData(2:nData+1);
  counter = counter + 1;
end % iWave
clear tData waveTrain waveTrainFFT tDiffData tDataFFT;

fprintf(1, '%s\n', 'Done.');

%% Widths, Depths
% Wavelet widths are based on the total area of the hat function that is
% above the zero line.  The derivative searcher looks for crossing points,
% and the upper cap is (hatWidth + 1)./3 wide.
depthStep = mean(diff(inputDepth));
waveletDerivativeWidths = (waveletWidth + 1)./3.*depthStep;

%% Find the different regions in the entire transform
% We look for regions in the transform that have all of one sign, since the
% wavelet blocker essentially divides layers based on their second
% derivative.  A positive region implies that the borehole curve is going
% through a positive deflection, whilst a negative sign in the data says
% the curve is tending to a negative deflection.  We look for these regions
% and base the layers upon this.
fprintf(1, '%s', 'Searching for layers and layer boundaries...');
diffDataSign = (sign(diffData));
[layerLabelMatrix, numberOfLayers] = boreilabel(diffDataSign, 4);
fprintf(1, '%s\n', 'Done.');
fprintf(1, '%s\n', ['Number of layers: ' num2str(numberOfLayers)]);


fprintf(1, '%s', 'Sorting layers based on importance (deflection)...');

waitCounter = numberOfLayers*3;
waitHandle = waitbar(0, 'Reconciling layers...');

%% Eliminate layers that don't get back to the LHS.
% Sometimes there is a layer that doesn't make it back to the LHS of
% the transform.  These orphans need to be placed in the surrounding layer.
nLayer = numberOfLayers;
existingLayers = [];
regionCounter = 1;
countSoFar = 0;
for iLayer = 1:nLayer
  countSoFar = countSoFar +1;
  waitbar(countSoFar./waitCounter, waitHandle);
  % First, we look for the layer existing in the first column of the
  % labelled wavelet transform matrix.
  isInLayer = find(layerLabelMatrix(:,1) == iLayer, 1, 'first');
  
  % If the layer does not exist in the fist column of the labelled first
  % column matrix, then we look for the region to the immediate left of the
  % lowest column index.  We then change the existing layer so that it is equal 
  % to the index of the layer to the left.  When we do this, we must
  % update the number of layers, the layer index matrix, and the layer area
  % matrix.
  if isempty(isInLayer)
%     iLayer
    findLayerIndex = find(layerLabelMatrix == iLayer);
    [r, c] = ind2sub(size(layerLabelMatrix), findLayerIndex);
    [minC, minCIndex] = min(c(:));
    newLayer = layerLabelMatrix(r(minCIndex), minC - 1);
    layerLabelMatrix(findLayerIndex) = newLayer;
    numberOfLayers = numberOfLayers-1;
  % Otherwise, we continue on...
  else
    existingLayers(regionCounter, 1) = iLayer;
    regionCounter = regionCounter + 1;
  end
end % iLayer
nLayer = numberOfLayers;
waitCounter = waitCounter./3 + nLayer*2;

%% Rename layers
% Then we go back into the layerLabelMatrix and rename the layers so that
% they are contiguous 1:numberOfLayers.  Also find the indexes of the
% layers in the matrix.
% We go back through the modified matrix and calculate the mean of the
% wavelet transform region determined by the layer boundaries.  This is a
% way of evaluating the importance of the layer.
counter = 1;
if nargout >=2
  layerMeanMatrix = zeros(size(layerLabelMatrix));
end

for iLayer = existingLayers'
  countSoFar = countSoFar +1;
  waitbar(countSoFar./waitCounter, waitHandle, 'Searching for indexes...');
  findLayers{counter} = find(layerLabelMatrix == iLayer);
  layerMean(counter,1) = abs(mean(diffData(findLayers{counter})));
  if nargout >=2
    layerMeanMatrix(findLayers{counter}) = layerMean(counter,1);
  end
  counter = counter + 1;
end % iLayer

% layerMean = layerMean./max(layerMean);
layerMean = (layerMean-max(layerMean))./(max(layerMean)-min(layerMean))+1;
if nargout >=2
  layerMeanMatrix = (layerMeanMatrix - max(layerMeanMatrix(:)))./ ...
    (max(layerMeanMatrix(:)) - min(layerMeanMatrix(:))) + 1;
end 
[layerSort, layerImportance] = sort(layerMean, 'descend');

layerLabelMatixTemp = zeros(size(layerLabelMatrix));
for iLayer = 1:numel(existingLayers)
  countSoFar = countSoFar +1;
  waitbar(countSoFar./waitCounter, waitHandle, 'Filling the matrix...');
  layerLabelMatixTemp(findLayers{layerImportance(iLayer)}) = iLayer;
end
layerLabelMatrix = layerLabelMatixTemp;
clear layerLabelMatixTemp layerMean findLayers;


fprintf(1, '%s\n', 'Done.');
close(waitHandle);

%% Sort the layers
layerDepthIndexSorted = zeros(nLayer, 2);
for iLayer = 1:nLayer
  layerDepthIndexSorted(iLayer, 1) = ...
    find(layerLabelMatrix(:,1) == iLayer, 1, 'first');
  layerDepthIndexSorted(iLayer, 2) = ...
    find(layerLabelMatrix(:,1) == iLayer, 1, 'last');  
end % iRegion

layerDepth = inputDepth(layerDepthIndexSorted);
layerThickness = diff(layerDepth, 1, 2) + depthStep;

%% Put data into outputs
outputData = struct(...
  'depth', inputDepth, ...
  'data', inputData, ...
  'depthStep', depthStep, ...
  'nData', nData, ...
  'nWavelet', nWavelet, ...
  'nLayer', nLayer, ...
  'waveletTransform', diffData, ...
  'layerLabel', layerLabelMatrix, ...
  'layerMean', layerSort, ...
  'layerDepthIndex', layerDepthIndexSorted, ...
  'layerDepth', layerDepth, ...
  'layerThickness', layerThickness, ...
  'waveletWidth', waveletDerivativeWidths);
