
clear;
close all;

addpath('../StateReconstruction/');

load ../StateReconstruction/resultsAvrg/PromoterFeatures.mat;


[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);

durations = durations(durations>5);

TotMat = [];
CondMat = [];

features = {'pActivate', 'meanTrOutput', 'meanTauA', 'meanTauS'};

featureNames = {};
for k=1:length(durations)
    [Mat, promVec, conditionMat] = GetDataMatrix(dataB, '', [], [1], durations(k), [-1], features);
    
    TotMat = [TotMat, Mat];
    
    for i=1:length(features)
        featureNames{end+1} = [features{i} ':' num2str(durations(k)) 'min']; 
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

rowCount = 1;
PermMat = zeros(length(featureNames), length(featureNames));

for k=1:length(features)
    featureIdx = find(startsWith(featureNames, features{k}));
    
    for l=1:length(featureIdx)
       PermMat(rowCount, featureIdx(l)) = 1;
       featureNamesPerm{rowCount} = featureNames{featureIdx(l)};
       rowCount = rowCount + 1;
    end
end

TotMatPerm = (PermMat*TotMat')';

cg = clustergram(TotMatPerm, 'Cluster', 'Col', 'Standardize', 'col');
set(cg, 'ColumnLabels', featureNamesPerm);
set(cg, 'RowLabels', promNames);

