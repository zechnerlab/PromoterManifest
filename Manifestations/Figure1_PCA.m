
clear;
close all;

addpath('../Data/Msn2');
addpath('../Data/');
addpath('../StateReconstruction');

load ../StateReconstruction/results/PromoterFeatures.mat;


[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

durations = durations(durations>5);

%durations = 40;

TotMat = [];
CondMat = [];

features = {'pActivate', 'maxTr', 'meanTauA', 'meanTauS', 'timeToMaxTr'};
featureNames = {};

for k=1:length(durations)
    [Mat, promVec, conditionMat] = GetDataMatrix(dataB, '', [], [1], durations(k), [-1], features);
    
    TotMat = [TotMat, Mat];
    
    for i=1:length(features)
        featureNames{end+1} = [features{i} ':' num2str(durations(k))]; 
    end
    
end

exclPromoters = {'pSIP18_mut6', 'pSIP18_mut21'};
for k=1:length(exclPromoters)
    delIdx = strcmp(promVec, exclPromoters{k});
TotMat = TotMat(~delIdx, :);
promVec = promVec(~delIdx, :);
conditionMat = conditionMat(~delIdx, :);
end


validIdx = sum(isnan(TotMat), 2)==0;
TotMat = TotMat(validIdx, :);
promVec = promVec(validIdx);
conditionMat = conditionMat(validIdx);

% group normalization
% for k=1:length(features)
%     featureIdx = startsWith(featureNames, features{k});
%     vals = TotMat(:, featureIdx);
%     %vals = zscore(vals);
%     meanV = mean(vals(:));
%     stdV = std(vals(:));
%     vals = (vals - meanV) / stdV;
%     TotMat(:, featureIdx) = vals;
% end

TotMat = zscore(TotMat, [], 1);
[loadings, scores] = pca(TotMat, 'Centered', false);

colVec = GetDefaultColors();%distinguishable_colors(length(promoters));

coloredPromoters = {'DDR2', 'DCS2', 'HXK1', 'SIP18', 'TKL2'};
%coloredPromoters = {'HXK1', 'TKL2', 'DDR2'};
for k=1:size(TotMat, 1)
    subplot(2,2,1);
    p = plot(scores(k, 1), scores(k, 2), 'o'); hold on;

    
    
    colIdx = find(strcmp(promoters, promVec{k}));
    size = 1.5*find(conditionMat(k, 1) == concentrations);
    
    if (sum(strcmp(promVec{k}, coloredPromoters))>0)
        dotColor = colVec(colIdx, :);
    else
       dotColor = [0.9, 0.9, 0.9]; 
    end
    
    set(p, 'MarkerEdgeColor', dotColor, 'MarkerFaceColor', dotColor);
    set(p, 'MarkerSize', 1 + size);
    
    subplot(2,2,2);
    p = plot(scores(k, 1), scores(k, 3), 'o'); hold on;
    
    
    %colIdx = find(strcmp(promoters, promVec{k}));
    %dotColor = colVec(colIdx, :);
    
    set(p, 'MarkerEdgeColor', dotColor, 'MarkerFaceColor', dotColor);
    set(p, 'MarkerSize', 1 + size);
    
    subplot(2,2,3);
    bar(loadings(:, 1));
    
    subplot(2,2,4);
    bar(loadings(:, 2));
    
    promStr{k} = promVec{k};
end

subplot(2, 2, 1);
%legend(promStr(1:length(promoters)));
xlabel('PC 1');
ylabel('PC 2');
%xlim([-8 8])
%ylim([-3 3]);
plot([0, 0], ylim, 'k--');
plot(xlim, [0, 0], 'k--');

subplot(2, 2, 2);
xlabel('PC 2');
ylabel('PC 3');
xlim([-15 15])
ylim([-10 10]);
plot([0, 0], ylim, 'k--');
plot(xlim, [0, 0], 'k--');


