
clear all;
close all;

addpath('../StateReconstruction/');

load ../StateReconstruction/results/PromoterFeatures.mat;


[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

durations = durations(durations>5);

%durations = 50;

TotMat = [];
CondMat = [];

%features = {'pActivate', 'maxTr', 'meanTrOutput', 'meanTauA', 'meanTauS', 'timeToMaxTr'};
features = {'pActivate', 'meanTrOutput', 'meanTauA', 'meanTauS'};

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

for i=1:length(promVec)
   MSN2Induction = (find(concentrations == conditionMat(i)))*100/4;
   promNames{i} = [promVec{i} ' (' num2str(MSN2Induction) '% Msn2)'];
end

% group normalization
for k=1:length(features)
    featureIdx = startsWith(featureNames, features{k});
    vals = TotMat(:, featureIdx);
    vals = zscore(vals);
    TotMat(:, featureIdx) = vals;
end

%TotMat = zscore(TotMat, [], 1);
%[loadings, scores] = pca(TotMat, 'Centered', false);

cg = clustergram(TotMat, 'Cluster', 'All', 'Standardize', 'none');
set(cg, 'ColumnLabels', featureNames);
set(cg, 'RowLabels', promNames);

