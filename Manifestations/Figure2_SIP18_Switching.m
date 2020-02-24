clear all;
close all;

addpath('results/');
addpath('../Data/MSN2');
addpath('../Common/');
addpath('../StateReconstruction/');

load ../StateReconstruction/results/PromoterFeatures.mat;


[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

durations = durations(durations>5);

durations = [10, 20, 30, 40, 50];

colVec = GetDefaultColors();
features = dataB.featureNames;
numFeatures = length(features);

featureIdx = nchoosek(1:numFeatures, 2);
featureIdx = unique(featureIdx, 'rows');



numPlotsPerFig = 5;
rowIdx = 0;
figWidth= 3000;
figHeight = 1000;


legendStr = {};

selPromoters = {'RTN2', 'SIP18', 'pSIP18_mut6', 'pSIP18_mut21'};

semFactor = 2.58;

for u=1:length(selPromoters)
    l = find(strcmp(promoters, selPromoters{u}));
    
    feature1 = 'TotalSwitches';
    sem1 = 'TotalSwitchesStd';

    p = zeros(length(concentrations), 1);
    
    concentrationsP = [25, 50, 75, 100];
    
    for k=1:length(durations)
        
        [data1, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], feature1);
        [semData1, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], sem1);
        
        data1 = data1;%./durations(k);
        semData1 = semData1;%./durations(k);
        
        figure(1);
        subplot(length(durations), length(selPromoters), (k-1)*length(selPromoters) + u);
        bar(concentrationsP, data1); hold on;
        errorbar(concentrationsP, data1, semFactor*semData1, '.');
        
        Dat(:, k) = data1;
    end
    %figure(2);
    %subplot(1, length(selPromoters), u);
    %heatmap(Dat, 'CellLabelColor', 'none', 'GridVisible', 'off', 'ColorbarVisible', 'on');
    %title(selPromoters{u});
    
    figure(3);
    subplot(1, length(selPromoters), u);
    plot(Dat');
end

