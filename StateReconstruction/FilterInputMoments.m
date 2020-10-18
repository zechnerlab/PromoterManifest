function [Mu, Sigma, PosteriorTime, Likelihood, PosteriorStats] = FilterInputMoments(model, config, t, Z, targetData, M0, measureIdx)

model.TranscriptionRateSequence.inputTimes = t(2:end);
model.TranscriptionRateSequence.inputLevels = Z(1:end);

if (nargin<5)
    M0 = model.InputODEInfos.infos.DefaultInitialConditions;
    means = [model.InitialMeans; model.muMeas0; model.H(1)];
    vars = [model.InitialVars; model.varMeas0; model.H(2)];

    M0 = SetInitialInputMoments(model.InputODEInfos.infos, M0, means, vars);
end

if (nargin<6)
   measureIdx = 1:length(targetData.MeasurementTime); 
end


numBins = model.NumBins;
numMoments = length(M0);

tPathY = [targetData.MeasurementTime(measureIdx)];
yPath = targetData.Measurement(measureIdx);


currT = tPathY(1);

PosteriorStats = []; %zeros(numMoments, (numBins - 1) * length(tPathY) - numBins);
PosteriorTime = [];%zeros(1, (numBins - 1) * length(tPathY) - numBins);

ODEOpts.RelTol = 0.01;

dT = config.ODEStepSize;

transcrSequence.t = t;
transcrSequence.Value = [Z];

for k=2:length(tPathY)
    nextT = tPathY(k);
    
    numBins = max(3, ceil((nextT-currT)/dT));
    
    time = linspace(currT, nextT, numBins);
    
    model.Update = 0;
    
    output = SolveGeneModel(time, model, config, transcrSequence, M0);
    
    if (output.valid==0)
        Mu = 0;
        Sigma = 0;
        Likelihood = -inf;
       return;
    end
    
    xPrior = output.M;
    time = output.t;
    if (numBins == 2)
        xFiltered = [xPrior(:, 1), xPrior(:, end)];
    else
        xFiltered = xPrior;
    end
    


    try
        model.Update = 1;
        model.Measurement = yPath(k);
        [posteriorStat, LTmp] = model.InputODEInfos.odeHandle(0, xFiltered(:, end), model);
        xFiltered(:, end) = xFiltered(:, end) + posteriorStat;
        L(:, k-1) = LTmp;
        
        M0 = xFiltered(:, end);
        idx = (k-2)*(numBins-1)+1:(k-1)*(numBins-1);
        PosteriorStats = [PosteriorStats, xFiltered(:, 1:end-1)];
        PosteriorTime = [PosteriorTime, time(1:end-1)];
        
        currT = nextT;
        
        if (model.StatusOutput)
            fprintf('Processed measurement %d...\n', measureIdx(k-1));
        end        
    catch msg
        Mu = 0;
        Sigma = 0;
        Likelihood = -inf;
        return;
    end
end



PosteriorTime(end+1) = time(end);
PosteriorStats(:, end+1) = xFiltered(:, end);

% compute prior probabilities for likelihood computation

Likelihood = L;

Mu = PosteriorStats(1:3, :);
Sigma = PosteriorStats(4:9, :);

end