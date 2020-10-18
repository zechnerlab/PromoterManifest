function [MAPParams, MMSEParams, chain] = RunMCMCMoments(model, config, Data, startParams, M, sigmaProp, validIdx)



numParams = length(startParams);

if (nargin<6)
    sigmaProp = 0.08;
end

chain = zeros(numParams, M);
chain(:, 1) = startParams;
bestL = -inf;

w = getCurrentTask;
if (~isempty(w))
    workerId = w.ID;
else
    workerId = 0;
end

for i=1:M
    
    currParams = chain(:, i);
    
    forwardProb = 0;
    backwardProb = 0;
    
    propParams = currParams;
    if (i>1)
        %propose new parameter set using lognormal distribution and
        %calculate forward and backward probability (needed for
        %Metropolis-Hastings).
        for u=1:numParams
            k = u;
            propParams(k) = lognrnd(log(currParams(k)), sigmaProp);
            forwardProb = forwardProb + log(lognpdf(propParams(k), log(currParams(k)), sigmaProp));
            backwardProb = backwardProb + log(lognpdf(currParams(k), log(propParams(k)), sigmaProp));
        end
        
    end
    
    %Evaluate the objective function (i.e., the Posterior density). This
    %function solves the moment equations and evaluates the likelihood
    %function and priors. LNew is the log-posterior density evaluated at
    %the proposed parameter set.
    [LNew, Stats, Time, currModel] = ObjectiveMoments(propParams,  model, config, Data, validIdx);
    
    %Calculate Metropolis-Hastings acceptance probability.
    if (i>1)    
        if (~isnan(LNew))
            a = min(1, exp(LNew - LOld + backwardProb - forwardProb));
        else
            a=0;
        end
    else
        a = 1;
    end
    
    %Accept with probability a.
    if (rand<a)
        
        LOld = LNew;
        chain(:, i+1) = propParams;
        if (LOld>bestL)
            MAPResults.Stats = Stats;
            MAPResults.Time = Time;
            MAPResults.Model = currModel;
            MAPResults.BestParams = propParams;
            bestL = LOld;
        end
    else
        chain(1:numParams, i+1) = chain(1:numParams, i);
    end
    
    
    %Store current Posterior density value.
    LVec(i) = LOld;
    
    %This part prints / plots the current status of the inference. Since
    %this step is slow, it is executed only every 200 iterations.
    if (mod(i, 200)==0)
        fprintf('MCMC worker (%d): it=%d, L=%f:', workerId, i, LOld);
        
        fprintf('cT = %f, Z(1) = %f, Z(2) = %f, Z(3)= %f, sigma=%f\n', MAPResults.Model.c(3), MAPResults.Model.Z(1), MAPResults.Model.Z(2), MAPResults.Model.Z(3), MAPResults.Model.MeasurementSigma);
        
        
        figure(37);
        subplot(2,4,1);
        semilogy(chain(:, 1:i)');
        
        subplot(2,4,2);
        plot(LVec(1:i));
        
        for l=1:length(MAPResults.Stats)
            Stats = MAPResults.Stats{l};
            Time = MAPResults.Time;

            mRNAMean = Stats(4, :);
            proteinMean = Stats(5, :);
            proteinVar = Stats(25, :) - proteinMean.^2;
            mRNAVar = Stats(22, :) - mRNAMean.^2;

            subplot(2,4,3);
            plot(Time, Stats([3], :)); hold on;

            subplot(2,4,4);
            plot(Time, Data.Conditions{l}.Moments(1, :), 'o', Time, proteinMean);  hold on;
            
            subplot(2,4,4);
            plot(Time, Data.Conditions{l}.Moments(1, :) + sqrt(Data.Conditions{l}.Moments(2, :)), '-', Time, Data.Conditions{l}.Moments(1, :) - sqrt(Data.Conditions{l}.Moments(2, :)), '-');
            plot(Time, proteinMean + sqrt(proteinVar), Time, proteinMean - sqrt(proteinVar));  hold on;
        
        end
        subplot(2,4,3);
        hold off;
        
        subplot(2,4,4);
        hold off;
        
        subplot(2,4,5);
        hold off;
        
        subplot(2,4,6);
        bar(MAPResults.BestParams);
        
        subplot(2,4,7);
        plot(Time, mRNAMean, Time, mRNAMean - sqrt(mRNAVar), Time, mRNAMean + sqrt(mRNAVar));
        
      
        drawnow;
    end
    
end

%Store and return inference results.
[~, mapIdx] = max(LVec);
MAPParams = chain(:, mapIdx);
%Burnin of 3000 samples is used. Not relevant for the MAP estimates though
%since those are just determined by finding the parameter set which
%achieved the highest posterior density. Relevant only for MMSE estimates
%(which are not used in this study).
burnIn = 3000;
MMSEParams = mean(chain(:, burnIn:end), 2);
