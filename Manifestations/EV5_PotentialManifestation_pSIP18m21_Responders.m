clear;
close all;

addpath('../StateReconstruction');
addpath('../Common/');
addpath('../Data/Msn2');

resultsName = 'resultsAvrg';

pathStr = ['../StateReconstruction/' resultsName '/PromoterFeatures.mat'];
load(pathStr);

[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

colVec = GetDefaultColors();


features = dataB.featureNames;
numFeatures = length(features);

featureIdx = nchoosek(1:numFeatures, 2);
featureIdx = unique(featureIdx, 'rows');


legendStr = {};

promName = 'pSIP18_mut21';

l = find(strcmp(promoters, promName));

feature1 = 'meanTauA';
feature2 = 'meanTrOutputR';
feature3 = 'meanTrOutput';

semFactor = 2;

p = zeros(length(concentrations), 1);
durations = durations(durations > 5);
for k=1:length(concentrations)
    
    
    [data1, promVec, conditionMat, data1sem] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature1);
    [data2, promVec, conditionMat, data2sem] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature2);
    [data3, promVec, conditionMat, data3sem] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature3);
   
    data1 = data1/60;
    data1sem = data1sem/60;

    figure(2);
    
    subplot(2, length(concentrations), k);
    
    p = bar(durations, data1, 0.7); hold on;
    
    %binofit(numResponders, numCells);
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(durations, data1, data1sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    xlim([0, 60]);
    title(num2str(concentrations(k)));
    ylabel('Time active');
    xlabel('Duration');
    
    subplot(2, length(concentrations), k + length(concentrations));
    

    p = errorbar(durations, data2, data2sem*semFactor, 'ro-'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    xlim([0, 60]);
    ylabel('Transcr output');
    xlabel('Duration');
    
    p = errorbar(durations, data3, data3sem*semFactor, 'bo-');
    set(p, 'MarkerFaceColor', 'b');
    set(p, 'MarkerEdgeColor', 'b');
    
    
end

