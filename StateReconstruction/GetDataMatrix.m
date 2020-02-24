function [Mat, promVec, conditionMat] = GetDataMatrix(dataB, promoter, conc, numPulses, duration, interval, features)

selIdx = ones(size(dataB.Data, 1), 1);

if (strcmp(features, ''))
    featureIdx = 1:length(dataB.featureNames);
else
    if (iscell(features))
    for k=1:length(features)
        featureIdx(k) = find(strcmp(dataB.featureNames, features{k}));
    end
    elseif (ischar(features))
       featureIdx = find(strcmp(dataB.featureNames, features));
    end
    
end


if (~strcmp(promoter, ''))
    selIdx = and(selIdx, strcmp(dataB.Promoter, promoter)');
end
if (~isempty(conc))
    selIdx = and(selIdx, dataB.Condition(:, 1) == conc);
end
if (~isempty(numPulses))
    selIdx = and(selIdx, dataB.Condition(:, 2) == numPulses);
end
if (~isempty(duration))
    selIdx = and(selIdx, dataB.Condition(:, 3) == duration);
end
if (~isempty(interval))
    selIdx = and(selIdx, dataB.Condition(:, 4) == interval);
else
    selIdx = and(selIdx, dataB.Condition(:, 4) ~= -1);  
end


Mat = dataB.Data(selIdx, featureIdx);
promVec = dataB.Promoter(selIdx)';
conditionMat = dataB.Condition(selIdx, :);
end