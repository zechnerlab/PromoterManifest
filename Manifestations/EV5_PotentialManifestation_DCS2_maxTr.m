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

promName = 'DCS2';
l = find(strcmp(promoters, promName));


feature1 = 'meanTauA';
feature2 = 'maxTr';

semFactor = 2;

p = zeros(length(concentrations), 1);

durations = durations(durations>5);
for k=1:length(concentrations)
    
    
    [data1, promVec, conditionMat, data1sem] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature1);
    [data2, promVec, conditionMat, data2sem] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature2);

    numCells  = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], 'numCells');
    numResponders = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], 'numResponders');
    
    data1 = data1/60;
    data1sem = data1sem/60;
    
    
    [datR] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], 'pActivate'); 
    
    figure(2);
    
    subplot(2, length(concentrations), k);
    
    p = bar(durations, data1, 0.7); hold on;
    
    %binofit(numResponders, numCells);
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(durations, data1, data1sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    xlim([0, 60]);
    ylim([0, 50]);
    title(num2str(concentrations(k)));
    ylabel('Time active (min)');
    xlabel('Duration');
    
    subplot(2, length(concentrations), k + length(concentrations));
    
    p = bar(durations, data2, 0.7); hold on;
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(durations, data2, data2sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    %ylim([0, 0.15]);
    xlim([0, 60]);
    ylim([0, 0.15]);
    ylabel('Max Transcr.');
    xlabel('Duration');
    
    legendStr{k} = num2str(concentrations(k));
end

