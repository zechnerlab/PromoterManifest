clear all;
close all;

addpath('results/');
addpath('../Data/MSN2');

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'TKL2'};

counter = 0;
maxCond = 8;


legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [30];

conditions = GenerateConditions();

colVec = GetDefaultColors();
colVec = flipud(colVec);

selPromoters = {'DCS2', 'DDR2', 'TKL2', 'SIP18', 'pSIP18_mut21'};
selPromoters = {'pSIP18_mut21'};

for i=1:length(selPromoters)
    
    k = find(strcmp(promoters, selPromoters{i}));
    
    cGrid = linspace(0, 20, 20);
    for u=20:-1:1

        concIdx = find(conditions{u}.Concentration == cVec);
        plIdx = floor(conditions{u}.PulseParameters(2)/10);
        pulseLength = conditions{u}.PulseParameters(2);
        if sum(pulseLength == plVec)==0
            continue;
        end
        if (sum(pulseLength == [])>0)
            continue;
        end
        
        results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
        
        responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
        %idx = find(responders);%find(results.Promoters{k}.Valid);
        idx = find(results.Promoters{k}.Valid);

        results.Promoters{k}.Model.Z
        
        [~, maxZIdx] = min(results.Promoters{k}.Model.Z(2:3));
        maxZIdx = maxZIdx + 1;
        meanTauTot{concIdx, plIdx} = sum(results.Promoters{k}.MeanTauState(idx, [2, 3]), 2);

        transcrOutput{concIdx, plIdx} = results.Promoters{k}.MeanTrOutput(idx);
   
    end

    
    numConc = size(meanTauTot, 1);
    numPL = size(meanTauTot, 2);
    
    for u=1:numConc
        
        meanTauTotVec = [];
        meanTranscrOutputVec = [];
        for l=1:numPL
            meanTauTotVec = [meanTauTotVec; meanTauTot{u, l}];
            meanTranscrOutputVec = [meanTranscrOutputVec; transcrOutput{u, l}];
        end
        
        %CorrVec(concIdx) = corr(meanTauA/60, maxTr);
        if (sum(results.Promoters{k}.Valid))
            m = fitlm(meanTauTotVec/60, meanTranscrOutputVec);
            CorrVecTauA(u) = m.Rsquared.Ordinary;
            SlopeVecTauA(u) = m.Coefficients.Estimate(2);
            
            %m = fitlm(meanTau2/60, meanTau3/60);
            %CorrVecTau23(concIdx) = m.Rsquared.Ordinary;
            %SlopeVecTau23(concIdx) = m.Coefficients.Estimate(2);
        else
            %CorrVecTau23(concIdx) = nan;
            %SlopeVecTau23(concIdx) = nan;
            
            CorrVecTauA(u) = nan;
            SlopeVecTauA(u) = nan;
        end
        
        rdIdx = randperm(length(meanTauTotVec));
        
        %rdIdx = rdIdx(1:80);
        figure(3);
        subplot(3, length(selPromoters), i);
        p = scatter(meanTauTotVec(rdIdx)/60, meanTranscrOutputVec(rdIdx), 'o'); hold on;
        pVec(u) = p;
        xlabel('Time active');
        ylabel('Transcr. output');
        set(p, 'SizeData', 15);
        %set(p, 'MarkerFaceColor', colVec(u, :));
        %set(p, 'MarkerEdgeColor', colVec(u, :));
        set(p, 'MarkerFaceColor', get(p, 'MarkerEdgecolor'));
        set(p, 'MarkerEdgeColor', 'none');
        set(p, 'MarkerFaceAlpha', 1);
        %set(p, 'MarkerEdgeAlpha', 0.2);
        title(selPromoters{i});
        box on;
        
        figure(3);
        subplot(3, length(selPromoters), length(selPromoters) + i);
        bar(CorrVecTauA);
        ylabel('R Squared');
        
        subplot(3, length(selPromoters), 2*length(selPromoters) + i);
        bar(SlopeVecTauA);
        ylabel('Slope');
        
        drawnow;
    end
    
    for u=numConc:-1:1
       uistack(pVec(u), 'top'); 
    end
end
    
    
