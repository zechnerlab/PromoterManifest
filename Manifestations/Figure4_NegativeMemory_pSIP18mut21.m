
clear;
close all;

addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');

%specify results to average over.
resultsNames{1} = 'results_1';
resultsNames{2} = 'results_2';
resultsNames{3} = 'results_3';
resultsNames{4} = 'results_4';
resultsNames{5} = 'results_5';

conditionsAll = GenerateConditions();

semFactor = 2;

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;

%Specify index of promoter to show (1=DCS2, 9=SIP18 D6);
kIdx = [1, 9];

%Only plot the repeated pulse conditions with different intervals
cIdx = 20+[3, 7, 8, 9, 10];

colVec = 'rgbkcmy';

for i=1:length(kIdx)
    k = kIdx(i);
    for j=1:length(cIdx)
        
        u = cIdx(j);
        condIdx = j;
        
        for l=1:length(resultsNames)
            resultsName = resultsNames{l};
            results = load(['../StateReconstruction/' resultsName '/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');

            numCells = size(results.Promoters{k}.P2, 1);
            
            trRateMat(l, :) = mean(results.Promoters{k}.MeanTr(logical(results.Promoters{k}.Valid), :), 1);
        end
        
        meanTr = mean(trRateMat, 1);
        q(1, :) = meanTr - semFactor*std(trRateMat, [], 1)/sqrt(length(resultsNames));
        q(2, :) = 2*semFactor*std(trRateMat, [], 1)/sqrt(length(resultsNames));
        
        t = results.Promoters{k}.tGrid;
        
        figure(1);
        subplot(length(cIdx), 2, (condIdx-1)*length(kIdx) + i);
        p = area(t/60, q'); hold on;
        set(p(1), 'FaceColor', 'none');
        set(p(1), 'EdgeColor', 'none');
        set(p(2), 'FaceAlpha', 0.2, 'FaceColor', colVec(1), 'EdgeColor', colVec(1));
        plot(t/60, meanTr, 'Color', colVec(1));
        xlabel('Time in min');
        ylabel('\lambda(t)');
        xlim([0, 100]);
        title([promoters{k} '(' num2str(conditions{u}.Name) ')'], 'interpreter', 'none');
        drawnow;
        
        if (i==1)
            figure(2);
            subplot(length(cIdx), 1, j);
            plot(results.Promoters{k}.YFP{1}.InputParams.inputTimes/60, results.Promoters{k}.YFP{1}.InputParams.inputLevels(1:end-1));
            xlim([0, 100]);
            xlabel('Time in min');
            ylabel('Msn2 level');
        end
        
        drawnow;
        
    end
end

