clear;
close all;

addpath('../StateReconstruction');
addpath('../Common/');
addpath('../Data/Msn2');

resultsName = 'resultsAvrg';

pathStr = ['../StateReconstruction/' resultsName '/PromoterFeatures.mat'];
load(pathStr);


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

[~, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

colVec = GetDefaultColors();

selPromoters = {'DCS2', 'DDR2', 'TKL2'};

feature1 = 'RSquared';
feature2 = 'beta';


semFactor = 2;
plotIndividualRepeats = 0;

p = zeros(length(concentrations), 1);

durations = durations(durations>5);

%specify which 
concVecPercent = [25 50 75 100]; %concentration in %
concVec = [100, 275, 690, 3000]; %concentration in nM
duration = 50;

for l=1:length(selPromoters)

    [data1, promVec, conditionMat, data1sem, runs1] = GetDataMatrix(dataB, selPromoters{l}, [], [1], duration, [-1], feature1);
    [data2, promVec, conditionMat, data2sem, runs2] = GetDataMatrix(dataB, selPromoters{l}, [], [1], duration, [-1], feature2);
    
    plotIndividualRepeats = and(plotIndividualRepeats, ~isempty(runs1));

    if (plotIndividualRepeats)
        runs1 = reshape(runs1, size(runs1, 1), size(runs1, 3));
        runs2 = reshape(runs2, size(runs2, 1), size(runs2, 3));
    end
    
    validIdx = 1:length(data1);
    
    subplot(2, length(selPromoters), l);
    
    p = bar(concVecPercent, data1, 0.7); hold on;
    
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(concVecPercent, data1, data1sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    ylim([0, 1]);
    xlim([0, 125]);
    ylabel('R squared');
    title(selPromoters{l});
    
    if (plotIndividualRepeats)
        plot(concVecPercent, runs1, '.');
    end
    
    
    
    subplot(2, length(selPromoters), l + length(selPromoters));
    
    p = bar(concVecPercent, data2, 0.7); hold on;
    set(p, 'FaceColor', [0.6, 0.6, 0.8]);
    p = errorbar(concVecPercent, data2, data2sem*semFactor, 'ro'); hold on;
    set(p, 'MarkerFaceColor', 'r');
    set(p, 'MarkerEdgeColor', 'r');
    ylabel('k');
    %ylim([0, 0.6]);
    xlim([0, 125]);
end

drawnow;

%% Plot one of the runs to show the individual single-cell reconstructions

conditions = GenerateConditions();

conditions = GetConditions(conditions, duration, concVec);

colVec = GetDefaultColors();
colVec = flipud(colVec);

resultsName = 'results_1'; %specify which run to show;

figure;
for i=1:length(selPromoters)
    
    k = find(strcmp(promoters, selPromoters{i}));
    
    cGrid = linspace(0, 20, 20);
    for u=1:length(conditions)

        concIdx = find(conditions{u}.Concentration == concVec);
        plIdx = find(conditions{u}.PulseParameters(2) == duration);

        pulseLength = conditions{u}.PulseParameters(2);

        results = load(['../StateReconstruction/' resultsName '/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');

        validIdx = results.Promoters{k}.Valid;
        idx = find(validIdx);

        [~, maxZIdx] = min(results.Promoters{k}.Model.Z(2:3));
        maxZIdx = maxZIdx + 1;
        meanTau = sum(results.Promoters{k}.MeanTauState(idx, [2, 3]), 2);
        transcrO = results.Promoters{k}.MeanTrOutput(idx);
        vIdx = meanTau<inf*60;
        meanTauTot{concIdx, plIdx} = meanTau(vIdx);
        transcrOutput{concIdx, plIdx} = transcrO(vIdx);

        subplot(length(duration), length(selPromoters), (plIdx - 1)*length(selPromoters) + i);
        p = scatter(meanTauTot{concIdx, plIdx}/60, transcrOutput{concIdx, plIdx}, 'o'); hold on;
        pVec(u) = p;
        xlabel('Time transcribing');
        ylabel('Transcr. output');
        set(p, 'SizeData', 15);
        set(p, 'MarkerFaceColor', get(p, 'MarkerEdgecolor'));
        set(p, 'MarkerEdgeColor', 'none');
        set(p, 'MarkerFaceAlpha', 1);
        title(selPromoters{i});
        box on;
        drawnow;
    end

    
   for u=length(conditions):-1:1
      uistack(pVec(u), 'top'); 
   end
end


