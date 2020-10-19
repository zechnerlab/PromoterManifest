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

legendStr = {};

promName = 'RTN2';
l = find(strcmp(promoters, promName));

feature1 = 'cv2TrOutput';
feature2 = 'meanTrOutput';


semFactor = 2;


p = zeros(length(concentrations), 1);

durations = [10, 50];
concVec = [25 50 75 100];

for k=1:length(durations)
    
    
    [data1, promVec, conditionMat, data1sem] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], feature1);
    [data2, promVec, conditionMat, data2sem] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], feature2);

    validIdx = 1:length(data1);

    figure(2);
    
    subplot(2, length(durations), k);
    
    p = bar(concVec, data1, 0.7); hold on;

    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(concVec, data1, data1sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    ylabel('Noise transcr. output');
    xlabel('Msn2 conc');
    
    subplot(2, length(durations), k + length(durations));
    
    p = bar(concVec, data2, 0.7); hold on;
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(concVec, data2, data2sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    ylabel('Average transcr. output');
    xlabel('Msn2 conc');

end

