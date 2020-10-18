function [PosteriorStats, PosteriorTime, Likelihood] = RunFilter(targetData, options)


cVec = options.c;

M0 = options.M0;
numBins = options.NumBins;
numMoments = length(M0);

tPathY = targetData.MeasurementTime;
yPath = targetData.Measurement;

currT = 0;
PosteriorStats = zeros(numMoments, (numBins - 1) * length(tPathY) - numBins);
PosteriorTime = zeros(1, (numBins - 1) * length(tPathY) - numBins);

ODEOpts.RelTol = 0.001;
ODEOpts.AbsTol = 0.001;

for k=2:length(tPathY)
    nextT = tPathY(k);
    
    time = linspace(currT, nextT, numBins);
    
    options.Update = 0;
    options.InputParams = targetData.InputParams;
    
    funHandleDop = @(t, y) options.odeInfos.odeHandle(t, y, options);
    
    
    [~, xPrior] = ode15s(funHandleDop, time, M0, ODEOpts);
    
    
    %[xPrior] = ode2(funHandleDop, time, M0);
    
    
    if (numBins == 2)
        xFiltered = [xPrior(1, :); xPrior(end, :)];
    else
        xFiltered = xPrior;
    end
    
    
    try
        options.Update = 1;
        options.Measurement = yPath(k);
        [posteriorStat, ~, LTmp] = options.odeInfos.odeHandle(0, xFiltered(end, :), options);
        xFiltered(end, :) = xFiltered(end, :) + posteriorStat;
        L(:, k-1) = LTmp;
        LTest(k-1) = sum(xFiltered(end,1:2).*LTmp);
        
        M0 = xFiltered(end, :)';
        idx = (k-2)*(numBins-1)+1:(k-1)*(numBins-1);
        PosteriorStats(:, idx) = xFiltered(1:end-1, :)';
        PosteriorTime(idx) = time(1:end-1);
        
        currT = nextT;
        
        if (options.StatusOutput)
            fprintf('Processed measurement %d...\n', k-2);
        end
        
    catch msg
        
        Likelihood = -inf;
        return;
    end
end

PosteriorTime(end+1) = time(end);
PosteriorStats(:, end+1) = xFiltered(end, :)';

% compute prior probabilities for likelihood computation

Likelihood = L;

end