clear all;
close all;

addpath('results/');

load ../StateReconstruction/results/PromoterFeatures.mat;


[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

colVec = GetDefaultColors();

features = dataB.featureNames;
numFeatures = length(features);

featureIdx = nchoosek(1:numFeatures, 2);
featureIdx = unique(featureIdx, 'rows');



numPlotsPerFig = 5;
rowIdx = 0;
figWidth= 3000;
figHeight = 1000;


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();
conditions = GetConditions(conditions, plVec, cVec);

colVec = 'rgbkcmy';

k=find(strcmp(promoters, 'SIP18'));
numRows = length(cVec);
numCols = length(plVec);

for u=1:length(conditions)
    
    results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
    %for k=1:length(results.Promoters)
        t = results.Promoters{k}.tGrid;
        
        concIdx = find(conditions{u}.Concentration == cVec);
        plIdx = find(conditions{u}.PulseParameters(2) == plVec);
        
        numCells = size(results.Promoters{k}.P2, 1);
        
        responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
        meanTr = mean(results.Promoters{k}.MeanTr(responders, :), 1);
        
        numValidCells = sum(results.Promoters{k}.Valid);
        
        r = mean(results.Promoters{k}.MeanRSwitches);
        Q = [0    r(3) r(5)
            r(1) 0    r(6)
            r(2) r(4) 0   ];
        
        
        figure(k);
        subplot(numRows, numCols, u);
        %HeatMap(Q, 'CellLabelColor', 'none', 'GridVisible', 'off', 'ColorbarVisible', 'off');
        p = image(Q);
        title(conditions{u}.Name, 'Interpreter', 'none');
        
        
        figure(99);
        p = plot(r(1)/(plIdx*10), r(6)/(plIdx*10), 'o'); hold on;
        set(p, 'MarkerFaceColor', colVec(plIdx), 'MarkerEdgeColor', colVec(plIdx));
        set(p, 'MarkerSize', 2+concIdx);
        drawnow;
end





