function [L, StatsOut, TimeOut, model] = ObjectiveMoments(commonParams, model, config, targetData, validIdx)

if sum(commonParams<0)>0
    L = -inf;
    StatsOut = 0;
    TimeOut = 0;
    model = 0;
    return;
end

model = SetModelParameters(model, config, commonParams);

if sum(model.P0<0)>0
    L = -inf;
    fprintf('Invalid initial condition\n');
    StatsOut = 0;
    TimeOut = 0;
    model = 0;
    return;
end

numConditions = length(targetData.Conditions);

TimeOut = targetData.Time;

%This loop is relevant only if multiple conditions are fitted at once.
%Here, only one condition is fitted at a time.
for k=1:numConditions
    
    %Solve moment equations using the current parameter set (stored in
    %'model').
    [Stats] = RunModel(model, targetData.Time, targetData.Conditions{k}.InputParams);
    
    mRNAMean = Stats(4, :);
    proteinMean = Stats(5, :);
    proteinVar = Stats(25, :) - proteinMean.^2;
    mRNAVar = Stats(22, :) - mRNAMean.^2;
    
    %Set the predicted moments and the moments obtained from data as well
    %as their uncertainties.
    targetMoments.Mean = targetData.Conditions{k}.Moments(1, :);
    targetMoments.Variance = targetData.Conditions{k}.Moments(2, :);
    targetUncertainties.Mean = targetData.Conditions{k}.Uncertainties(1, :);
    targetUncertainties.Variance = targetData.Conditions{k}.Uncertainties(2, :);
    targetUncertainties.MeanVariance = targetData.Conditions{k}.Uncertainties(3, :);
    
    PredictedMoments.Mean = proteinMean;
    PredictedMoments.Variance = proteinVar + model.MeasurementSigma^2;
    
    %Calculate Loglikelihood
    LVec(k) = EvaluateLogLikelihood(targetMoments, targetUncertainties, PredictedMoments, validIdx);
    
    StatsOut{k} = Stats;
end

%Add the logarithm of the prior distribution. The Prior is calculated in
%'EvaluatePrior'.
L = sum(LVec) + EvaluatePrior(model, config);

%Some sanity checks.
if (~isreal(L))
    L = -inf;
end

if (isinf(L))
    L = -inf;
end



end