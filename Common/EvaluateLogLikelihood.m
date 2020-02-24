function [L] = EvaluateLogLikelihood(TargetMoments, TargetUncertainties, PredictedMoments, idx)

    MomentsT(1, :) = TargetMoments.Mean;
    MomentsT(2, :) = TargetMoments.Variance;

    MomentsP(1, :) = PredictedMoments.Mean;
    MomentsP(2, :) = PredictedMoments.Variance;
    
    Uncertainties(1, :) = TargetUncertainties.Mean;
    Uncertainties(2, :) = TargetUncertainties.Variance;
    Uncertainties(3, :) = TargetUncertainties.MeanVariance;
    
    %Uncertainties = max(1000, Uncertainties);% = eps;
    
    for i=1:length(idx)
        x = MomentsT(:, idx(i));
        y = MomentsP(:, idx(i));
        Sigma = [Uncertainties(1, idx(i)), Uncertainties(3, idx(i))
                 Uncertainties(3, idx(i)), Uncertainties(2, idx(i))];
        LPart(i) = -1/2 * (x - y)' * (Sigma \ (x - y));
    
    end
    
    L = sum(LPart);

end