clear all;
close all;


addpath('../Data/MSN2');
addpath('../Common');

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'TKL2'};

counter = 0;
maxCond = 8;


legendStr = {};
cVec = [100, 3000];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();



colVec = GetDefaultColors();
%colVec = flipud(colVec);

selPromoter = 'pSIP18_mut21';

k = find(strcmp(promoters, selPromoter));


cGrid = linspace(0, 20, 20);
for i=1:length(plVec)
    for j=1:length(cVec)
        
        condition = GetConditions(conditions, plVec(i), cVec(j));
        
        results = load(['../StateReconstruction/results/StateReconstruction_' condition{1}.Name '.mat'], 'Promoters');
        
        responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
        idxResponders = find(responders);%find(results.Promoters{k}.Valid);
        idxAll = find(results.Promoters{k}.Valid);
        idxNonResponders = setdiff(idxAll, idxResponders);
        
        [~, maxZIdx] = min(results.Promoters{k}.Model.Z(2:3));
        maxZIdx = maxZIdx + 1;
        meanTauTotAll{j, i} = sum(results.Promoters{k}.MeanTauState(idxAll, [2, 3]), 2);
        transcrOutputAll{j, i} = results.Promoters{k}.MeanTrOutput(idxAll);
        
        meanTauTotR{j, i} = sum(results.Promoters{k}.MeanTauState(idxResponders, [2, 3]), 2);
        transcrOutputR{j, i} = results.Promoters{k}.MeanTrOutput(idxResponders);
        
        meanTauTotNR{j, i} = sum(results.Promoters{k}.MeanTauState(idxNonResponders, [2, 3]), 2);
        transcrOutputNR{j, i} = results.Promoters{k}.MeanTrOutput(idxNonResponders);
        
        Responders{j, i} = length(idxResponders)/length(idxAll);
        
%         
%         meanTauTotAll{j, i} = sum(results.Promoters{k}.MeanTauA(idxAll), 2);
%         transcrOutputAll{j, i} = results.Promoters{k}.MeanTrOutput(idxAll);
%         
%         meanTauTotR{j, i} = sum(results.Promoters{k}.MeanTauA(idxResponders), 2);
%         transcrOutputR{j, i} = results.Promoters{k}.MeanTrOutput(idxResponders);
%         
%         meanTauTotNR{j, i} = sum(results.Promoters{k}.MeanTauA(idxNonResponders), 2);
%         transcrOutputNR{j, i} = results.Promoters{k}.MeanTrOutput(idxNonResponders);
%         
        
        mean(meanTauTotR{j, i})

        subplot(1,5,i);
%         p = scatter(meanTauTotNR{j, i}/60, transcrOutputNR{j, i}, 'o'); hold on;
%         xlabel('Time active');
%         ylabel('Transcr. output');
%         set(p, 'SizeData', 15);
%         set(p, 'MarkerFaceColor', colVec(i, :));
%         set(p, 'MarkerEdgeColor', 'none');
%         set(p, 'MarkerFaceAlpha', 0.05);
%         set(p, 'MarkerEdgeAlpha', 0.05);
%         title(selPromoter);
%         box on;
        
       
        p = scatter(meanTauTotAll{j, i}/60, transcrOutputAll{j, i}, 'o'); hold on;
        xlabel('Time active');
        ylabel('Transcr. output');
        set(p, 'SizeData', 15);
        set(p, 'MarkerFaceColor', colVec(j, :));
        set(p, 'MarkerEdgeColor', 'none');
        set(p, 'MarkerFaceAlpha', 1);
        set(p, 'MarkerEdgeAlpha', 1);
        title(selPromoter);
        box on;
        p = plot(mean(meanTauTotAll{j, i})/60, mean(transcrOutputAll{j, i}), 'v');
        set(p, 'MarkerSize', 10);
        set(p, 'MarkerFaceColor', colVec(j, :));
        set(p, 'MarkerEdgeColor', 'k');
        
        p = plot(mean(meanTauTotR{j, i})/60, mean(transcrOutputR{j, i}), 'd');
        set(p, 'MarkerSize', 10);
        set(p, 'MarkerFaceColor', colVec(j, :));
        set(p, 'MarkerEdgeColor', 'k');
        xlim([0, 60]);
        ylim([0, 800]);
        
%         subplot(1, 2, 2);
%         [h, b] = hist(meanTauTotR{j, i}/60);
%         h = h/(sum(h)*(b(2)-b(1)));
%         
%         area(b, h, 'FaceColor', colVec(i, :), 'FaceAlpha', 0.2); hold on;
        
        drawnow;
        
    end
end



