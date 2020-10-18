clear;
close all;

addpath('../Data/MSN2');
addpath('../Common/');

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];
colVec = GetDefaultColors();

conditions = GenerateConditions();


selPromoters = {'DCS2', 'TKL2', 'SIP18', 'DDR2', 'pSIP18_mut6', 'pSIP18_mut21'};

count = 1;

actThreshold = 0;

resultsName = 'resultsAvrg';

pathStr = ['../StateReconstruction/' resultsName '/PromoterFeatures.mat'];
load(pathStr);

[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

feature1 = 'meanTrOutput';
feature2 = 'cv2TrOutput';

for u=1:30
    
    concIdx = find(conditions{u}.Concentration == cVec);
    plIdx = floor(conditions{u}.PulseParameters(2)/10)+1;
    pulseLength = conditions{u}.PulseParameters(2);
    
    if (length(conditions{u}.PulseParameters) > 2)
       pulseInterval =  conditions{u}.PulseParameters(3);
    else
       pulseInterval = -1; 
    end
    
    for l = 1:length(selPromoters)
    [data1, promVec, conditionMat, data1sem] = GetDataMatrix(dataB, selPromoters{l}, ...
        conditions{u}.Concentration, conditions{u}.PulseParameters(1), pulseLength, pulseInterval, feature1);
    [data2, promVec, conditionMat, data2sem] = GetDataMatrix(dataB, selPromoters{l}, ...
        conditions{u}.Concentration, conditions{u}.PulseParameters(1), pulseLength, pulseInterval, feature2);

    if (sum(pulseLength == [])>0)
        continue;
    end
    
    if (plIdx>1)
        markerType = 'o';
    else
        markerType = 'v';
    end
    
    %results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
        figure(1);
        subplot(3, length(selPromoters), l);
        p = loglog(data1, data2, markerType); hold on;
        set(p, 'Color', colVec(plIdx, :));
        set(p, 'MarkerFaceColor', colVec(plIdx, :));
        set(p, 'MarkerSize', 2 + concIdx*2);
        xlabel('Average');
        ylabel('CV2');
        title(selPromoters{l});
        %xlim([1, 1000]);
        %ylim([0.1, 10]);
        
        
        subplot(3, length(selPromoters), length(selPromoters) + l);
        p = loglog(data1, data2, markerType); hold on;
        set(p, 'Color', colVec(concIdx, :));
        set(p, 'MarkerFaceColor', colVec(concIdx, :));
        set(p, 'MarkerSize', 2 + plIdx*2);
        xlabel('Average');
        ylabel('CV2');
        title(selPromoters{l});
        %xlim([1, 1000]);

%         subplot(3, length(selPromoters), 2*length(selPromoters) + i);
%         p = loglog(mean(proteinOutput), (std(proteinOutput)/mean(proteinOutput))^2, markerType); hold on;
%         set(p, 'Color', colVec(concIdx, :));
%         set(p, 'MarkerFaceColor', colVec(concIdx, :));
%         set(p, 'MarkerSize', 2 + plIdx*2);
%         xlabel('Average');
%         ylabel('CV2');
%         title(selPromoters{i});
%         %xlim([1, 1000]);
        drawnow;
    end
end

% 
% for i=1:length(selPromoters)    
%     subplot(2, length(selPromoters), i);
%     xLims = xlim;
%     muVec = linspace(xLims(1), xLims(2), 100);
%     plot(muVec, 1./sqrt(muVec));
%     
%     
%     subplot(2, length(selPromoters), length(selPromoters) + i);
%     xLims = xlim;
%     muVec = linspace(xLims(1), xLims(2), 100);
%     plot(muVec, 1./sqrt(muVec));
%     drawnow;
% end
% 
