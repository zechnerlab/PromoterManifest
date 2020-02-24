clear;
close all;

warning off;

addpath('../Common');
addpath('../Common/Statistics/');
addpath('../Common/StochChemKin/');
addpath('../Common/Models/');
addpath('../Common/ODEs/dopri');
addpath('../Common/ODEs/');
addpath('../Common/');
addpath('../Common/SpecialFunctions/');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');



conditionsAll = GenerateConditions();
%conditions = conditions(21:26);


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'DCS2'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;

colVec = 'rgbkcmy';


figWidth= 300;
figHeight = 3000;
%fig = figure('rend','painters','pos',[1000 1000 figWidth figHeight]);

kIdx = [1, 9];%find(strcmp(promoters, {'DCS2', 'pSIP18_mut21'}));

cIdx = 20+[3, 7, 8, 9, 10];

for i=1:length(kIdx)
    k = kIdx(i);
    for l=1:length(cIdx)
        
        u = cIdx(l);
        
        
        results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
        
        concIdx = find(conditions{u}.Concentration == cVec);
        numPulseIdx = floor(conditions{u}.PulseParameters(1));
        condIdx = l;
        
        
        
        numCells = sum(results.Promoters{k}.Valid);% size(results.Promoters{k}.P2, 1);
        
        responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
        valids = logical(results.Promoters{k}.Valid);
        meanTr = mean(results.Promoters{k}.MeanTr(valids, :));
        stdTr = std(results.Promoters{k}.MeanTr(valids, :));
        
        f = 2.58;
        q(1, :) = meanTr - f*std(results.Promoters{k}.MeanTr(valids, :))./sqrt(sum(valids));
        q(2, :) = 2*f*std(results.Promoters{k}.MeanTr(valids, :))./sqrt(sum(valids));
        
        t = results.Promoters{k}.tGrid;
        
        figure(1);
        subplot(length(cIdx), 2, (condIdx-1)*length(kIdx) + i);
        p = area(t/60, q'); hold on;
        set(p(1), 'FaceColor', 'none');
        set(p(1), 'EdgeColor', 'none');
        set(p(2), 'FaceAlpha', 0.2, 'FaceColor', colVec(1), 'EdgeColor', colVec(1));
        stairs(t/60, meanTr, 'Color', colVec(1));
        xlabel('Time in min');
        ylabel('\lambda(t)');
        xlim([0, 100]);
        ylim([-0.01, 0.1])
        title([promoters{k} '(' num2str(conditions{u}.Name) ')']);
        drawnow;
        
        if (i==1)
            figure(2);
            subplot(length(cIdx), 1, l);
            plot(results.Promoters{k}.YFP{1}.InputParams.inputTimes/60, results.Promoters{k}.YFP{1}.InputParams.inputLevels(1:end-1));
            xlim([0, 100]);
            xlabel('Time in min');
            ylabel('Msn2 level');
        end
        
        %scan for max transcription in each pulse
        window = zeros(size(meanTr));
        
        numPulses = conditions{u}.PulseParameters(1);
        pulseLength = conditions{u}.PulseParameters(2);
        pulseInterval = conditions{u}.PulseParameters(3);
        maxTrPerPulse = [];
        
        activationPulse = ceil(results.Promoters{k}.MeanTauS(responders) / ((pulseLength+pulseInterval)*60));
        cellsActivePerPulse = [];
        numValidCells = sum(responders);%numCells;
        for j=1:numPulses
            minT = (j-1)*(pulseLength+pulseInterval)*60;
            maxT = minT + (pulseLength+pulseInterval)*60;
            window = find(and(t>=minT, t<maxT));
            
            [maxTrPerPulse(j), maxIdx] = max(meanTr(window));
            %maxTrPerPulse(l) = mean(max(results.Promoters{k}.MeanTr(valids, window)));
            stdOfMaxTrPerPulse(j) = stdTr(window(maxIdx));
            
            cellsActivePerPulse(j) = sum(activationPulse==j)/numCells;%/numValidCells;
            [~, bounds] = binofit(sum(activationPulse==j), numCells);
            cellsActiveLowerP(j) = bounds(1);
            cellsActiveUpperP(j) = bounds(2);
            %numValidCells = numValidCells - sum(activationPulse==l);
        end
        
        figure(20+k);
        errorbar(maxTrPerPulse, f*stdOfMaxTrPerPulse/sqrt(sum(valids)), 'o-'); hold on;
        xlabel('Pulse');
        ylabel('Max Tr');
        legendStr{l} = num2str(pulseInterval);
        legend(legendStr);
        
        
        figure(40);
        subplot(length(cIdx), 2, (condIdx-1)*length(kIdx) + i);
        bar(cellsActivePerPulse);
        %errorbar(1:4, cellsActivePerPulse, cellsActivePerPulse-cellsActiveLowerP, cellsActiveUpperP - cellsActivePerPulse);
        
        
        
        xlabel('Pulse')
        ylabel('% cells switched on');
        title(['p=' num2str(numValidCells/numCells)]); 
       
        
        drawnow;
        
    end
end


