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
    
    for k=1:numConditions
        options = model;
        [Stats] = RunModel(model, targetData.Time, targetData.Conditions{k}.InputParams);
        
        %remove time zero
        Stats = Stats(:, 1:end);
        
        numMoments = length(model.odeInfos.infos.MomentSystem{1}.dM);
        numStates = length(model.Z);
        
        baseIdx = (0:numStates-1)*numMoments + numStates;
        
        mRNAMean = sum(Stats(baseIdx+1, :));
        proteinMean = sum(Stats(baseIdx+2, :));
        proteinVar = sum(Stats(baseIdx+7, :)) - proteinMean.^2;
        mRNAVar = sum(Stats(baseIdx+4, :)) - mRNAMean.^2;
        
        targetMoments.Mean = targetData.Conditions{k}.Moments(1, :);
        targetMoments.Variance = targetData.Conditions{k}.Moments(2, :);
        targetUncertainties.Mean = targetData.Conditions{k}.Uncertainties(1, :);
        targetUncertainties.Variance = targetData.Conditions{k}.Uncertainties(2, :);
        targetUncertainties.MeanVariance = targetData.Conditions{k}.Uncertainties(3, :);

        PredictedMoments.Mean = proteinMean;
        PredictedMoments.Variance = proteinVar + model.MeasurementSigma^2;
        %PredictedMoments.Variance = sqrt((proteinVar+proteinMean.^2)*(1+model.MeasurementSigma^2) - proteinMean.^2);
        %PredictedMoments.Variance = sqrt((proteinVar+proteinMean.^2)*exp(2*model.MeasurementSigma^2) - proteinMean.^2*exp(model.MeasurementSigma^2));
        
        LVec(k) = EvaluateLogLikelihood(targetMoments, targetUncertainties, PredictedMoments, validIdx);
        
        StatsOut{k} = Stats;
    end
    L = sum(LVec) + EvaluatePrior(model, config);
    
    if (~isreal(L))
       L = -inf; 
    end
    
    if (isinf(L))
       L = -inf; 
    end
    
    
   
end