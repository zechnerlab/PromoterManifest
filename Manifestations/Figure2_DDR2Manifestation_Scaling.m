clear all;
close all;

addpath('results/');
addpath('../StateReconstruction');
addpath('../Common/');
addpath('../Data/Msn2');

resultsName = 'resultsAvrg';
pathStr = ['../StateReconstruction/' resultsName '/PromoterFeatures.mat'];
load(pathStr);


[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

colVec = GetDefaultColors();

plotIndividualRepeats = 0;


promName = 'DDR2';
l = find(strcmp(promoters, promName));

feature1 = 'meanTrOutputR';
feature2 = 'meanTauS';
feature3 = 'meanTauA';
feature4 = 'maxTr';
feature5 = 'pActivate';


semFactor = 2;

p = zeros(length(concentrations), 1);
durations = durations(durations >5);

for k=1:length(concentrations)
    
    
    [data1, promVec, conditionMat, data1sem, runs1] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature1);
    [data2, promVec, conditionMat, data2sem, runs2] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature2);
    [data3, promVec, conditionMat, data3sem, runs3] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature3);
    [data4, promVec, conditionMat, data4sem, runs4] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature4);
    [data5, promVec, conditionMat, data5sem, runs5] = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], feature5);
    numCells  = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], 'numCells');
    numResponders = GetDataMatrix(dataB, promoters{l}, concentrations(k), [1], [], [-1], 'numResponders');
    
    data1 = data1;
    data2 = data2/60;
    data3 = data3/60;
    data4 = data4;
    data5 = data5
    
    data2sem = data2sem/60;
    data3sem = data3sem/60;

    plotIndividualRepeats = and(plotIndividualRepeats, ~isempty(runs1));

    if (plotIndividualRepeats)
        runs1 = reshape(runs1, size(runs1, 1), size(runs1, 3));
        runs2 = reshape(runs2, size(runs2, 1), size(runs2, 3));
        runs3 = reshape(runs3, size(runs3, 1), size(runs3, 3));
        runs4 = reshape(runs4, size(runs4, 1), size(runs4, 3));
        runs5 = reshape(runs5, size(runs5, 1), size(runs5, 3));
    end
 
    validIdx = logical(~isnan(data1));
    
    figure(1);
    subplot(1, 2, 1);
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
    
    if (plotIndividualRepeats == 1)
        p = plot(runs1, runs2/60, '.');
        set(p, 'MarkerFaceColor', colVec(k, :), 'MarkerEdgeColor', colVec(k, :)); 
        set(p, 'MarkerSize', 2 + j*1);
    end

    subplot(1, 2, 2);
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
    
    if (plotIndividualRepeats == 1)
        p = plot(runs1, runs3/60, '.');
        set(p, 'MarkerFaceColor', colVec(k, :), 'MarkerEdgeColor', colVec(k, :)); 
    end
    
    figure(2);
    
    subplot(2, length(concentrations), k);
    
    p = bar(durations, data5, 0.7); hold on;
    
%    binofit(numResponders, numCells);
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    %p = errorbar(durations, data4, data4sem*semFactor, 'ro'); hold on;
    %set(p, 'MarkerFaceColor', 'r');
    %set(p, 'MarkerEdgeColor', 'r');
    ylim([0, 1]);
    xlim([0, 60]);
    
    subplot(2, length(concentrations), k + length(concentrations));
    
    p = bar(durations, data4, 0.7); hold on;
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(durations, data4, data4sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    ylim([0, 0.6]);
    xlim([0, 60]);
    
    legendStr{k} = num2str(concentrations(k));
end

