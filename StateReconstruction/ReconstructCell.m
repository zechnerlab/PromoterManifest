function particles = ReconstructCell(model, config, particles, cellStruct, filterLength)


%This function executes the hybrid sequential Monte Carlo algorithm.

M = length(particles);
numMeasurements = length(cellStruct.Measurement);
numRecursions = ceil(numMeasurements / (filterLength-1))-1;
measurementTimes = cellStruct.MeasurementTime;


LAll = ones(M, 1);

for l=1:numRecursions
    
    %determine the next update time intervals.
    measureIdx = ((l-1)*(filterLength-1) + 1):(min(l*(filterLength-1)+1, numMeasurements));
    nextT = measurementTimes(measureIdx(end));
    
    %some sanity checks.
    LAll(or(isnan(LAll), isinf(LAll))) = -inf;
    
    if (sum(abs(imag(LAll))>0)>0)
        LAll(abs(imag(LAll))>0) = -inf;
    end
    
    %normalize particle weights.
    W = exp(LAll - max(LAll));
    W = W/sum(W);
    if (sum(isnan(W))>0 || sum(W<0)>0)
        particles = {};
        return;
    end
    
    %re-draw particles according to their weights.
    idx = randsample(1:M, M, true, W);
    particles = particles(idx);
    
    %Forward simulate the promoter until the next time. This is done for
    %all particules simulatenously.
    [particles] = SimulatePromoter(model, particles, nextT);
    
    %For all particles, evaluate the marginal likelihood function to
    %determine their weight.
    for k=1:M

        %Use conditional moments associate to particle k as initial
        %condition.
        M0 = particles{k}.PosteriorStats(:, end);
        
        %Solve marginal likelihood
        [Mu, Sigma, PosteriorTime, Likelihood, PosteriorStats] = FilterInputMoments(particles{k}.Model, config, particles{k}.t, particles{k}.Z, cellStruct, M0, measureIdx);
        
        %Store the conditional moments
        particles{k}.PosteriorStats = [particles{k}.PosteriorStats, PosteriorStats(:, 2:end)];
        particles{k}.PosteriorTime = [particles{k}.PosteriorTime, PosteriorTime(2:end)];
        particles{k}.L = sum(Likelihood);
        LAll(k) = particles{k}.L;
        
    end
    
    %some sanity checks (again, just to be really safe).
    LAll(or(isnan(LAll), isinf(LAll))) = -inf;
    
    if (sum(abs(imag(LAll))>0)>0)
        LAll(abs(imag(LAll))>0) = -inf;
    end
    
    %Normalize particle weights.
    W = exp(LAll - max(LAll));
    W = W/sum(W);
    
    for k=1:M
        particles{k}.W = W(k);
    end
end



