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
features.cv2TrOutput = var(promoter.MeanTrOutput(v, :))/features.meanTrOutput^2;

%v = logical(promoter.Valid);
v = responders;
features.meanTrOutputR = mean(promoter.MeanTrOutput(v, :));
features.meanTrOutputRStd = sqrt(var(promoter.MeanTrOutput(v, :))/sum(v));
features.cv2TrOutputR = var(promoter.MeanTrOutput(v, :))/features.meanTrOutputR^2;

meanTau = sum(promoter.MeanTauState(logical(promoter.Valid), [2, 3]), 2);
m = fitlm(meanTau/60, promoter.MeanTrOutput(logical(promoter.Valid), :));
features.RSquared = m.Rsquared.Ordinary;
features.beta = m.Coefficients.Estimate(2);






end