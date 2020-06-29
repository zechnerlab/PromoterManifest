clear;
close all;

addpath('results/');
addpath('../Data/MSN2');
addpath('../Common/');

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'DCS2'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();

%colVec = 'rgbkcmy';
%col = [0.7, 1, 0.6];
%colVec = col.*[0.2; 0.4; 0.6; 0.8; 1];
colVec = GetDefaultColors();

figWidth= 1000;
figHeight = 1000;
%fig = figure('rend','painters','pos',[1000 1000 figWidth figHeight]);
figure;

%k=find(strcmp(promoters, 'pSIP18_mut21'));

selPromoters = {'RTN2'};
%selPromoters = {'SIP18'};
%selPromoters = promoters;
count = 1;

actThreshold = 0;

for u=1:20
    
    
    
    concIdx = find(conditions{u}.Concentration == cVec);
    plIdx = floor(conditions{u}.PulseParameters(2)/10)+1;
    pulseLength = conditions{u}.PulseParameters(2);
    
    
    if (sum(pulseLength == [])>0)
        continue;
    end
    
    results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
    for i=1:length(selPromoters)
        k = find(strcmp(promoters, selPromoters{i}));
        responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
        transcrOutput = results.Promoters{k}.MeanTrOutput(find(results.Promoters{k}.Valid));
        
        
        Y = [];
        count = 1;
        for j=1:length(results.Promoters{k}.Valid)
            if results.Promoters{k}.Valid(j)==1
                %figure(99);
                %plot(max(results.Promoters{k}.YFP{j}.Measurement), results.Promoters{k}.MeanTrOutput(j), 'o'); hold on;
                Y(count) = max(results.Promoters{k}.YFP{j}.Measurement);
                count = count + 1;
            end
            
        end
        %hold off;
        proteinOutput = Y;
        
        
        if sum(results.Promoters{k}.Valid) < actThreshold
            continue;
        end
        
        
        if (plIdx>1)
            markerType = 'o';
            noiseResults{i}.cv2Mat(concIdx, plIdx-1) = (std(transcrOutput)/mean(transcrOutput))^2;
            noiseResults{i}.meanMat(concIdx, plIdx-1) = mean(transcrOutput);
            
        else
            markerType = 'v';
        end
        
        
        
        figure(1);
        subplot(3, length(selPromoters), i);
        p = loglog(mean(transcrOutput), (std(transcrOutput)/mean(transcrOutput))^2, markerType); hold on;
        set(p, 'Color', colVec(plIdx, :));
        set(p, 'MarkerFaceColor', colVec(plIdx, :));
        set(p, 'MarkerSize', 2 + concIdx*2);
        xlabel('Average');
        ylabel('CV2');
        title(selPromoters{i});
        %xlim([1, 1000]);
        %ylim([0.1, 10]);
        
        
        subplot(3, length(selPromoters), length(selPromoters) + i);
        p = loglog(mean(transcrOutput), (std(transcrOutput)/mean(transcrOutput))^2, markerType); hold on;
        set(p, 'Color', colVec(concIdx, :));
        set(p, 'MarkerFaceColor', colVec(concIdx, :));
        set(p, 'MarkerSize', 2 + plIdx*2);
        xlabel('Average');
        ylabel('CV2');
        title(selPromoters{i});
        %xlim([1, 1000]);
        
        subplot(3, length(selPromoters), 2*length(selPromoters) + i);
        p = loglog(mean(proteinOutput), (std(proteinOutput)/mean(proteinOutput))^2, markerType); hold on;
        set(p, 'Color', colVec(concIdx, :));
        set(p, 'MarkerFaceColor', colVec(concIdx, :));
        set(p, 'MarkerSize', 2 + plIdx*2);
        xlabel('Average');
        ylabel('CV2');
        title(selPromoters{i});
        %xlim([1, 1000]);
        drawnow;
    end
end

close all;
figure;
subplot(2,2,1);
bar(noiseResults{1}.cv2Mat(:, 1));

subplot(2, 2, 2);
bar(noiseResults{1}.cv2Mat(:, end));


subplot(2,2,3);
bar(noiseResults{1}.meanMat(:, 1));

subplot(2, 2, 4);
bar(noiseResults{1}.meanMat(:, end));

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
