clear all;
close all;

addpath('../StateReconstruction');
addpath('../Common/');
addpath('../Data/MSN2');
addpath('../Data/');

colVec = GetDefaultColors();
semFactor = 2;

%specify results to average over.
resultsNames{1} = 'results_1';
resultsNames{2} = 'results_2';
resultsNames{3} = 'results_3';
resultsNames{4} = 'results_4';
resultsNames{5} = 'results_5';

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
promName = 'RTN2';

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 3000];
plVec = [10, 30, 50];

conditions = GenerateConditions();
conditions = GetConditions(conditions, plVec, cVec);


k=find(strcmp(promoters, promName));


for u=1:length(conditions)
    
    concIdx = find(conditions{u}.Concentration == cVec);
    plIdx = find(plVec == conditions{u}.PulseParameters(2));
    
    validRun = zeros(1, length(resultsNames));

    for l=1:length(resultsNames)
        resultsName = resultsNames{l};
        results = load(['../StateReconstruction/' resultsName '/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
        t = results.Promoters{k}.tGrid;
        
        numCells = size(results.Promoters{k}.P2, 1);
        responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
        trRateMat(l, :) = mean(results.Promoters{k}.MeanTr(responders, :), 1);
        trStdMat(l, :) = std(results.Promoters{k}.MeanTr(responders, :), [], 1);
        
        validRun(l) = sum(responders)>0;
    end
    
    validRun = logical(validRun);
    meanTr = mean(trRateMat(validRun, :), 1);
    stdTr = mean(trStdMat(validRun, :), 1);
    
    q(1, :) = meanTr - semFactor*std(trRateMat(validRun, :), [], 1)/sqrt(sum(validRun));
    q(2, :) = 2*semFactor*std(trRateMat(validRun, :), [], 1)/sqrt(sum(validRun));
    
    figure(99);
    subplot(length(cVec), 1, concIdx);
    
    validTimeIdx = find(t<100*60); %select time points before 100min for better vizualiation.
    p = area(t(validTimeIdx)/60, q(:, validTimeIdx)'); hold on;
    set(p(1), 'FaceColor', 'none');
    set(p(1), 'EdgeColor', 'none');
    set(p(2), 'FaceAlpha', 0.1, 'FaceColor', colVec(plIdx, :), 'EdgeColor', colVec(plIdx, :));
    plot(t(validTimeIdx)/60, meanTr(validTimeIdx), 'Color', colVec(plIdx, :));
    xlabel('Time in min');
    ylabel('\lambda(t)');
    %xlim([0, 100]);
    title([promoters{k} '(' num2str(conditions{u}.Concentration) 'nM)']);

    drawnow;
    
    fprintf('Processed promoter %s, condition %s\n', promoters{k}, conditions{u}.Name);
    
end
