clear all;
close all;

addpath('results/');
addpath('../Data/MSN2');
addpath('../Common/');

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

selPromoters = {'RTN2', 'ALD3', 'SIP18', 'pSIP18_mut6', 'pSIP18_mut21'};

semFactor = 2.58;

for u=1:length(selPromoters)
    l = find(strcmp(promoters, selPromoters{u}));
    
    feature1 = 'MeanSwitches21';
    sem1 = 'TotalSwitchesStd';
    feature2 = 'MeanSwitches23';
    sem2 = 'MeanActivationsStd';
    
    
    p = zeros(length(concentrations), 1);
    
    concentrationsP = [25, 50, 75, 100];
    
    for k=1:length(durations)
        
        [data1, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], feature1);
        [semData1, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], sem1);
        
        [data2, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], feature2);
        [semData2, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, [], [1], durations(k), [-1], sem2);
        
        data1 = data1;%./durations(k);
        semData1 = semData1./durations(k);
        
        data2 = data2;%./durations(k);
        semData2 = semData2./durations(k);
        
        
        Dat1(:, k) = data1;
        Dat2(:, k) = data2;
    end
    
    subplot(1,length(selPromoters), u);
    plot(Dat1', Dat2', 'o-'); drawnow;
end

