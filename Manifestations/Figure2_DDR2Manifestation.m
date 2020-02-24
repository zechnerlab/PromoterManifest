clear all;
close all;

addpath('results/');
addpath('../StateReconstruction');
addpath('../Common/');
addpath('../Data/Msn2');

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

promName = 'DDR2';
l = find(strcmp(promoters, promName));

feature1 = 'meanTrOutput';
feature2 = 'meanTauS';
feature3 = 'meanTauA';

sem1 = 'meanTrOutputStd';
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


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 690, 3000];
plVec = [10, 30, 50];

conditions = GenerateConditions();
conditions = GetConditions(conditions, plVec, cVec);
%colVec = 'rgbkcmy';

k=find(strcmp(promoters, promName));

%k = 2;
for u=1:length(conditions)
    
    results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    t = results.Promoters{k}.tGrid;
    
    
    concIdx = find(conditions{u}.Concentration == cVec);
    plIdx = find(plVec == conditions{u}.PulseParameters(2));
    
    numCells = size(results.Promoters{k}.P2, 1);
    
    responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
    %responders = logical(results.Promoters{k}.Valid);
    meanTr = mean(results.Promoters{k}.MeanTr(responders, :), 1);
    
    
    trRec = results.Promoters{k}.TrReconstruction(responders);
    
    numValidCells = sum(results.Promoters{k}.Valid);
    
    
    q(1, :) = meanTr - semFactor*std(results.Promoters{k}.MeanTr(responders, :), 1)/sqrt(sum(responders));
    q(2, :) = 2*semFactor*std(results.Promoters{k}.MeanTr(responders, :), 1)/sqrt(sum(responders));
    
    %q = quantile(results.Promoters{k}.MeanTr(responders, :), [0.3, 0.7], 1);
    
    
    figure(99);
    subplot(length(cVec), 1+length(plVec), (concIdx-1)*(1+length(plVec)) + 1);
    p = area(t/60, q'); hold on;
    set(p(1), 'FaceColor', 'none');
    set(p(1), 'EdgeColor', 'none');
    set(p(2), 'FaceAlpha', 0.1, 'FaceColor', colVec(plIdx, :), 'EdgeColor', colVec(plIdx, :), 'EdgeAlpha', 0.2);
    plot(t/60, meanTr, 'Color', colVec(plIdx, :));
    xlabel('Time in min');
    ylabel('\lambda(t)');
    xlim([0, 100]);
    title([promoters{k} '(' num2str(conditions{u}.Concentration) 'nM)']);
    
    
    subplot(length(cVec), 1+length(plVec), (concIdx-1)*(1+length(plVec)) + 1 + plIdx);
    
    cellIdx = randperm(numCells);
    responderData = [];
    nonresponderData = [];
    for i=1:length(cellIdx)
        j = cellIdx(i);
        
        %meanTr = results.Promoters{k}.MeanTr(j, :);
        %t = results.Promoters{k}.tGrid;
        
        %stairs(t/60, meanTr);
        
        %subplot(plotsPerFigure, plotsPerFigure*2, 2*counter);
        if (results.Promoters{k}.Valid(j)==1)
            
            %g = plot(results.Promoters{k}.YFP{j}.MeasurementTime/60, results.Promoters{k}.YFP{j}.Measurement, '-'); hold on;
            
            if (responders(j)==1)
                %    set(g, 'Color', [1, 0, 0, 0.5]);
                responderData = [responderData; results.Promoters{k}.YFP{j}.Measurement];
            else
                %    set(g, 'Color', [0, 0, 0, 0.5]);
                nonresponderData = [nonresponderData; results.Promoters{k}.YFP{j}.Measurement];
            end
        end
    end
    
    
    
    if (size(nonresponderData, 1)>0)
        qnonResponders = quantile(nonresponderData, [0.2, 0.8], 1);
        qnonResponders(2, :) = qnonResponders(2, :) - qnonResponders(1, :);
        
        %         p = area(results.Promoters{k}.YFP{1}.MeasurementTime/60, qnonResponders'); hold on;
        %         set(p(1), 'FaceColor', 'none');
        %         set(p(1), 'EdgeColor', 'none');
        %         set(p(2), 'FaceAlpha', 0.8*(size(nonresponderData, 1)/numValidCells), 'FaceColor', 'k', 'EdgeColor', 'none');
        
        n = size(nonresponderData, 1);
        idx = 1:min(5, n);
        plot(results.Promoters{k}.YFP{1}.MeasurementTime/60, nonresponderData(idx, :), 'r'); hold on;
    end
    %xlim([0, 50]);
    
    if (size(responderData, 1)>0)
        qResponders = quantile(responderData, [0.1, 0.9], 1);
        qResponders(2, :) = qResponders(2, :) - qResponders(1, :);
        
        %         p = area(results.Promoters{k}.YFP{1}.MeasurementTime/60, qResponders'); hold on;
        %         set(p(1), 'FaceColor', 'none');
        %         set(p(1), 'EdgeColor', 'none');
        %         set(p(2), 'FaceAlpha', 0.1 + 0*0.8*(size(responderData, 1)/numValidCells), 'FaceColor', 'r', 'EdgeColor', 'none');
        
        n = size(responderData, 1);
        %idx = randperm(n);
        idx = 1:min(5, n);
        plot(results.Promoters{k}.YFP{1}.MeasurementTime/60, responderData(idx, :), 'r'); hold on;
        %ylim([0, 5000]);
    end
    
    
    
    
    
    drawnow;
    
    fprintf('Processed promoter %s, condition %s\n', promoters{k}, conditions{u}.Name);
    
end



%
%
