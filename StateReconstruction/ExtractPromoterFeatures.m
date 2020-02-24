function [features] = ExtractPromoterFeatures(promoter, validOnly, maxTimeWindow)

responders = and(promoter.PActivated>0.99, promoter.Valid);
responders = and(responders, promoter.MeanTauS<maxTimeWindow);

features.pActivate = sum(responders) / sum(promoter.Valid);

if (validOnly==0)
   responders = and(ones(size(responders)), promoter.Valid); 
end

features.meanTauS = mean(promoter.MeanTauS(responders));
features.varTauS = var(promoter.MeanTauS(responders));
features.meanTauSStd = sqrt(features.varTauS/sum(responders));

features.meanTauA = mean(promoter.MeanTauA(responders));
features.varTauA = var(promoter.MeanTauA(responders));
features.meanTauAStd = sqrt(features.varTauA/sum(responders));


MeanTr = promoter.MeanTr(responders, :);
numResponders = sum(responders);

features.numResponders = numResponders;
features.numCells = sum(promoter.Valid);

if (numResponders>0)
for j=1:1000
    rdIdx = randi(numResponders, numResponders, 1);
    [maxTr, maxIdx] = max(mean(MeanTr(rdIdx, :)));
    maxTrVal(j) = maxTr;
    maxTrTime(j) = promoter.tGrid(maxIdx);
end

features.timeToMaxTr = mean(maxTrTime);
features.timeToMaxTrStd = std(maxTrTime);
features.maxTr = mean(maxTrVal);
features.maxTrStd = std(maxTrVal);

else
    features.timeToMaxTr = nan;
    features.maxTr = nan;
    features.maxTrStd = nan;
    features.timeToMaxTrStd = nan;
end
v = logical(promoter.Valid);
%v = responders;
features.meanTrOutput = mean(promoter.MeanTrOutput(v, :));
features.meanTrOutputStd = sqrt(var(promoter.MeanTrOutput(v, :))/sum(v));
features.varTrOutput = var(promoter.MeanTrOutput(v, :));

%v = logical(promoter.Valid);
v = responders;
features.meanTrOutputR = mean(promoter.MeanTrOutput(v, :));
features.meanTrOutputRStd = sqrt(var(promoter.MeanTrOutput(v, :))/sum(v));
features.varTrOutputR = var(promoter.MeanTrOutput(v, :));

%numCells = min(size(promoter.MeanTauState, 1), length(responders));
%idx = 1:numCells;

%features.MeanSwitches12 = mean(promoter.MeanRSwitches(idx, 1));
%features.MeanSwitches21 = mean(promoter.MeanRSwitches(idx, 3));
%features.MeanSwitches23 = mean(promoter.MeanRSwitches(idx, 4));
%features.MeanSwitches32 = mean(promoter.MeanRSwitches(idx, 6));

%features.VarSwitches12 = var(promoter.MeanRSwitches(idx, 1));
%features.VarSwitches21 = var(promoter.MeanRSwitches(idx, 3));
%features.VarSwitches23 = var(promoter.MeanRSwitches(idx, 4));
%features.VarSwitches32 = var(promoter.MeanRSwitches(idx, 6));


%features.MeanTimeState2 = mean(promoter.MeanTauState(idx, 2));
%features.MeanTimeState3 = mean(promoter.MeanTauState(idx, 3));
%features.VarTimeState2 = var(promoter.MeanTauState(idx, 2));
%features.VarTimeState3 = var(promoter.MeanTauState(idx, 3));

%features.c12 = mean(promoter.MeanRSwitches(idx, 1)) ./ mean(promoter.MeanTauState(idx, 1));
%features.c21 = mean(promoter.MeanRSwitches(idx, 3)) ./ mean(promoter.MeanTauState(idx, 2));
%features.c23  = mean(promoter.MeanRSwitches(idx, 4)) ./ mean(promoter.MeanTauState(idx, 2));
%features.c32  = mean(promoter.MeanRSwitches(idx, 6)) ./ mean(promoter.MeanTauState(idx, 3));

%features.TotalSwitches = mean(sum(promoter.MeanRSwitches(idx, [3, 6]), 2));
%features.TotalSwitchesStd = sqrt(var(sum(promoter.MeanRSwitches(idx, [3, 6]), 2)) / numCells);

%features.MeanActivations = mean(promoter.MeanNumSwitches(idx));
%features.MeanActivationsStd = sqrt(var(promoter.MeanNumSwitches(idx))/numCells);


end