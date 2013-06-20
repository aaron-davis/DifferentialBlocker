% DavisAndChristensenExample
clear all;
clc;
close all;

load('boreData');

%% Block the conduuctivity and gamma data
[blockData.cond.blockedMatrix, mMatrix, wMatrix] = BoreholeWaveletBlocker(...
  blockData.cond.data(:,2), blockData.cond.data(:,1));
[blockData.gamma.blockedMatrix, mMatrix, wMatrix] = BoreholeWaveletBlocker(...
  blockData.gamma.data(:,2), blockData.gamma.data(:,1));

%% Display the blocking results
% Conductivity
figure(1);
clf;
fHandle = DisplayBoreholeWaveletBlockerTransform(...
  blockData.cond.blockedMatrix,0,0,1);

% Gamma
figure(2);
clf;
fHandle = DisplayBoreholeWaveletBlockerTransform(...
  blockData.gamma.blockedMatrix,0,0,2);
  
% Layer importance
figure(3);
set(gcf, 'renderer', 'zbuffer');
pHandle = pcolor(blockData.gamma.blockedMatrix.waveletWidth, ...
  blockData.gamma.blockedMatrix.depth, ...
  log10(mMatrix));
shading flat;
colormap(flipud(hot(32)));
colorbar;
xlabel('Operator width (m)');
ylabel('Depth (m)');
title('Layer importance');

%% Display diagnostics
DisplayBoreholeWaveletBlockerDiagnostics(blockData.gamma.blockedMatrix, 4);

nLayerVal = 3;
percentageVal = 25;
minThicknessVal = 1;
operatorThicknessVal = 5;

layerBlocknLayers = BoreholeWaveletBlockerLayerSelection(...
  blockData.gamma.blockedMatrix, nLayerVal, 1);
layerBlockLayerPercent = BoreholeWaveletBlockerLayerSelection(...
  blockData.gamma.blockedMatrix, percentageVal, 2);
layerBlockLayerMinThick = BoreholeWaveletBlockerLayerSelection(...
  blockData.gamma.blockedMatrix, minThicknessVal , 3);
layerBlockLayerWaveletWidth = BoreholeWaveletBlockerLayerSelection(...
  blockData.gamma.blockedMatrix, operatorThicknessVal);

figure(5);
clf;
plot(blockData.gamma.data(:,2), blockData.gamma.data(:,1), 'k');
set(gca, 'ydir', 'reverse');
hold on;
plot(layerBlocknLayers.plotLayerMean, layerBlocknLayers.plotLayerDepth, 'r');
plot(layerBlockLayerPercent.plotLayerMean, ...
  layerBlockLayerPercent.plotLayerDepth, 'g');
plot(layerBlockLayerMinThick.plotLayerMean, ...
  layerBlockLayerMinThick.plotLayerDepth, 'b--');
plot(layerBlockLayerWaveletWidth.plotLayerMean, ...
  layerBlockLayerWaveletWidth.plotLayerDepth, 'm--');

xlabel('Counts');
ylabel('Depth (m)');
title('Natural gamma')









