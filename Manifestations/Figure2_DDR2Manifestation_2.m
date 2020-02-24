clear all;
close all;

addpath('results/');
addpath('../StateReconstruction');
addpath('../Common/');

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

%f = figure('rend','painters','pos',[1000 1000 figWidth figHeight]);

legendStr = {};

promName = 'DCS2';
l = find(strcmp(promoters, promName));

feature1 = 'maxTr';
feature2 = 'meanTauS';
feature3 = 'meanTauA';

sem1 = 'maxTrStd';
sem2 = 'meanTauSStd';
sem3 = 'meanTauAStd';

semFactor = 2.58;
%semFactor = 1;

p = zeros(length(concentrations), 1);

for k=1:length(concentrations)
    
    
    [data1, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature1);
    [data2, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature2);
    [data3, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature3);
    
    data1 = data1;
    data2 = data2/60;
    data3 = data3/60;
    
    [data1sem, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], sem1);
    [data2sem, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], sem2);
    [data3sem, promVec, conditionMat] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], sem3);
    
    data1sem = data1sem;
    data2sem = data2sem/60;
    data3sem = data3sem/60;
    
    [data4] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], 'pActivate'); 
    %data3 = data2 + ;
    
    %only allow if there's more than 5% cells responding.
    validIdx = data4>0.05;
    if (sum(validIdx)==0)
       continue; 
    end
    
    
    subplot(1, 3, 1);
    p(k) = plot(data1(validIdx), data2(validIdx), '-'); hold on;
    set(p(k), 'Color', colVec(k, :));
    Plot2DErrorbar(data1(validIdx), data2(validIdx), data1sem(validIdx), data2sem(validIdx), [colVec(k, :), 0.8], semFactor);
    
    for j=1:length(data1)
        if (validIdx(j)==1)
        p2 = plot(data1(j), data2(j), 'o');
        set(p2, 'MarkerFaceColor', colVec(k, :), 'MarkerEdgeColor', colVec(k, :));
        set(p2, 'MarkerSize', 2 + j*1);
        end
    end
    
    title(promoters{l});
    xlabel(feature1);
    ylabel(feature2);
    
    subplot(1, 3, 2);
    p(k) = plot(data1(validIdx), data3(validIdx), '-'); hold on;
    set(p(k), 'Color', colVec(k, :));
    Plot2DErrorbar(data1(validIdx), data3(validIdx), data1sem(validIdx), data3sem(validIdx), [colVec(k, :), 0.8], semFactor);
   
    
    for j=1:length(data1)
        if (validIdx(j)==1)
        p2 = plot(data1(j), data3(j), 'o');
        set(p2, 'MarkerFaceColor', colVec(k, :), 'MarkerEdgeColor', colVec(k, :));
        set(p2, 'MarkerSize', 2 + j*1);
        end
    end
    
    title(promoters{l});
    xlabel(feature1);
    ylabel(feature3);
    
    subplot(1,3,3);
    plot(durations, data4, 'o-'); hold on;

    legendStr{k} = num2str(concentrations(k));
end
