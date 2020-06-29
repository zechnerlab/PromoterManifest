function [features] = ExtractPromoterFeatures(promoter, validOnly)

responders = and(promoter.PActivated>0.99, promoter.Valid);
features.pActivate = sum(responders) / sum(promoter.Valid);

if (validOnly==0)
   responders = and(ones(size(responders)), promoter.Valid); 
end

features.meanTauS = mean(promoter.MeanTauS(responders));
features.varTauS = var(promoter.MeanTauS(responders));

features.meanTauA = mean(promoter.MeanTauA(responders));
features.varTauA = var(promoter.MeanTauA(responders));


[features.maxTr, maxIdx] = max(mean(promoter.MeanTr(responders, :)));
features.timeToMaxTr = promoter.tGrid(maxIdx);

features.meanTrOutput = mean(promoter.MeanTrOutput(responders, :));

numCells = min(size(promoter.MeanTauState, 1), length(responders));
idx = 1:numCells;

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

features.c12 = mean(promoter.MeanRSwitches(idx, 1)) ./ mean(promoter.MeanTauState(idx, 1));
features.c21 = mean(promoter.MeanRSwitches(idx, 3)) ./ mean(promoter.MeanTauState(idx, 2));
features.c23  = mean(promoter.MeanRSwitches(idx, 4)) ./ mean(promoter.MeanTauState(idx, 2));
features.c32  = mean(promoter.MeanRSwitches(idx, 6)) ./ mean(promoter.MeanTauState(idx, 3));

features.TotalSwitches = mean(sum(promoter.MeanRSwitches(idx, [3, 6]), 2));

features.MeanActivations = mean(promoter.MeanNumSwitches(idx));



end