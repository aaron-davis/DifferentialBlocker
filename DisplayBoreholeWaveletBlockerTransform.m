function funHandle = DisplayBoreholeWaveletBlockerTransform(inputData, ...
  matrixSwitch, plotTraceSwitch, figHandle)
%DisplayBoreholeWaveletBlockerTransform. Display the output of the wavelet blocking transform.
%   funHandle = DisplayBoreholeWaveletBlockerTransform(inputData).
%   Outputs the handle of the figure created to display the results of the
%   wavelet blocking transform.  The figure also contains the output of the
%   contour function whereby contours are drawn at the lavel of contour =
%   0. The contours show outlines of regions made by the wavelet transform.
%
%   funHandle = DisplayBoreholeWaveletBlockerTransform(inputData, matrixSwitch).
%   Same as above, but with the addition of a matrixSwitch flag that allows
%   the output of the regions created from the label function.  Regions
%   derived by the label function are the true layers of the wavelet
%   blocker, since the regions created by the blocker may not necessarily
%   go back to the LHS axis.  Any layers that do not go back to LHS of the
%   output are not true layers and are therefore sawlloed into layers that
%   contain them.  There are presently 2 switches for matrixSwitch:
%     matrixSwitch = 0 (default if no 2nd argument is input)
%       Displays the result from the wavelet blocker
%     matrixSwitch = 1
%       Displays the result of the label function, layers are coloured
%       according to the lines colourmap.
%     matrixSwitch = 2
%       Displays the result of the label function, layers are coloured
%       according to importance.
%
%
%   Company: Commonwealth Scientific and Industrial Research Organisation
%   (CSIRO), Earth Science and Resource Engineering, 2013
%   Author: Aaron C Davis
%
%   This software is licenced under the Creative Commons Attribution
%   (CC-BY) 3.0 licence (http://creativecommons.org/licenses/by/3.0/)

%% Detect the matrixSwitch
if nargin < 3
  plotTraceSwitch = 1;
  figHandle = 999;
end

%% Size the figure
PPI = get(0, 'screenPixelsPerInch');
if plotTraceSwitch
  pWidth = 5;
  pHeight = 5;
else
  pHeight = 5;
  pWidth = 3.3125;
end
  height = pHeight*PPI;
  width = pWidth*PPI;

%% Plot the figure
funHandle = figure(figHandle);
clf;
set(funHandle, 'Position', [0.5*PPI 0.5*PPI width height]);
set(funHandle, 'PaperUnits', 'inches');
set(funHandle, 'PaperSize', [pWidth pHeight]);
set(funHandle, 'PaperPosition', [0 0 pWidth pHeight]);
set(funHandle, 'Renderer', 'zbuffer');
set(funHandle, 'Color', 'w');

if plotTraceSwitch
  aHandle(1) = subplot(1,2,1);
  aHandle(2) = subplot(1,2,2);
else
  aHandle(1) = axes;
end
  
set(gcf, 'CurrentAxes', aHandle(1));

% Decide on the output, based on matrixSwitch
switch matrixSwitch
  case 0
    fieldName{1} = 'waveletTransform';
    cMap = [
                         0                         0                       0.5
        0.0355555555555555        0.0355555555555555         0.533333333333333
        0.0755555555555555        0.0755555555555555         0.566666666666667
                      0.12                      0.12                       0.6
         0.168888888888889         0.168888888888889         0.633333333333333
         0.222222222222222         0.222222222222222         0.666666666666667
                      0.28                      0.28                       0.7
         0.342222222222222         0.342222222222222         0.733333333333333
         0.408888888888889         0.408888888888889         0.766666666666667
                      0.48                      0.48                       0.8
         0.555555555555555         0.555555555555555         0.833333333333333
         0.635555555555556         0.635555555555556         0.866666666666667
                      0.72                      0.72                       0.9
         0.808888888888889         0.808888888888889         0.933333333333333
         0.902222222222222         0.902222222222222         0.966666666666667
                         1                         1                         1
         0.966666666666667         0.902222222222222         0.902222222222222
         0.933333333333333         0.808888888888889         0.808888888888889
                       0.9                      0.72                      0.72
         0.866666666666667         0.635555555555556         0.635555555555556
         0.833333333333333         0.555555555555555         0.555555555555555
                       0.8                      0.48                      0.48
         0.766666666666667         0.408888888888889         0.408888888888889
         0.733333333333333         0.342222222222222         0.342222222222222
                       0.7                      0.28                      0.28
         0.666666666666667         0.222222222222222         0.222222222222222
         0.633333333333333         0.168888888888889         0.168888888888889
                       0.6                      0.12                      0.12
         0.566666666666667        0.0755555555555555        0.0755555555555555
         0.533333333333333        0.0355555555555555        0.0355555555555555
                       0.5                         0                         0];
    titleString = 'Differential transform of data';
    cAxis = [-1 1]*max(abs(inputData.(fieldName{1})(:)));
  case 1
    fieldName{1} = 'layerLabel';
    cMapTemplate = get(0, 'FactoryAxesColorOrder');
    cMapTemplate =[cMapTemplate; 0.5 0.5 0.5];
    cMap = repmat(cMapTemplate, floor(inputData.nLayer./size(cMapTemplate, 1)), 1);
    cMap = [cMap; cMapTemplate(1:mod(inputData.nLayer, size(cMapTemplate,1)),:)];
    titleString = 'Layer regions';
    cAxis = [1 size(cMap, 1)];
  case 2
    fieldName{1} = 'layerLabel';
    cMap = [
                         0          0.47843137254902         0.752941176470588
                         0         0.493611638203669          0.76280834914611
                         0         0.508791903858318         0.772675521821632
                         0         0.523972169512966         0.782542694497154
                         0         0.539152435167615         0.792409867172676
                         0         0.554332700822264         0.802277039848197
        0.0221378874130298         0.566350411132195         0.808981657179001
        0.0487033523086654         0.577735610373182          0.81505376344086
        0.0752688172043011         0.589120809614168          0.82112586970272
         0.101834282099937         0.600506008855155         0.827197975964579
         0.128399746995572         0.611891208096142         0.833270082226439
         0.177735610373182         0.627324478178368         0.840860215053764
         0.238456672991777         0.644781783681214          0.84920936116382
         0.299177735610373         0.662239089184061         0.857558507273877
         0.359898798228969         0.679696394686907         0.865907653383934
         0.420619860847565         0.697153700189753         0.874256799493991
         0.472232764073371         0.716888045540797         0.882985452245414
         0.514737507906388         0.738899430740038         0.892093611638204
         0.557242251739405         0.760910815939279         0.901201771030993
         0.599746995572423          0.78292220113852         0.910309930423783
         0.642251739405439         0.804933586337761         0.919418089816572
         0.683997469955724         0.825426944971537         0.926502213788741
         0.724225173940544         0.842884250474383         0.929538266919671
         0.764452877925364          0.86034155597723         0.932574320050601
         0.804680581910183         0.877798861480076         0.935610373181531
         0.844908285895003         0.895256166982922          0.93864642631246
         0.880834914611006          0.91157495256167         0.942188488298545
         0.895256166982922         0.922201138519924         0.948260594560405
         0.909677419354839         0.932827324478178         0.954332700822264
         0.924098671726755         0.943453510436433         0.960404807084124
         0.938519924098672         0.954079696394687         0.966476913345984
         0.952941176470588         0.964705882352941         0.972549019607843];      
%     cMapTemplate = get(0, 'FactoryAxesColorOrder');
%     cMapTemplate =[cMapTemplate; 0.5 0.5 0.5];
%     cMap = repmat(cMapTemplate, floor(inputData.nLayer./size(cMapTemplate, 1)), 1);
%     cMap = [cMap; cMapTemplate(1:mod(inputData.nLayer, size(cMapTemplate,1)),:)];
    titleString = 'Layer importance';
    cAxis = [1 inputData.nLayer];
end % switch

pcolor(inputData.waveletWidth, inputData.depth, inputData.(fieldName{1}));
shading flat;
set(funHandle, 'colormap', cMap);
caxis(cAxis);
hold on;
% Place contours at 0, if wavelet transfer is selected
if matrixSwitch == 0
  contour(inputData.waveletWidth, inputData.depth, inputData.waveletTransform, ...
    [0 0], 'k', 'lineWidth', 0.72);
end % if matrixSwitch

% Annotate
set(gca, 'yDir', 'reverse');
ylabel('Depth (m)');
xlabel('Operator width (m)');
title(titleString);

if plotTraceSwitch
% Plot the log trace
set(gcf, 'CurrentAxes', aHandle(2));
plot(inputData.data, inputData.depth, 'k');
yLim = ylim;
set(gca, 'yDir', 'reverse');
ylabel('Depth (m)');
xlabel('Data value'); 
title('Trace');

set(gcf, 'CurrentAxes', aHandle(1));
ylim(yLim);
end;