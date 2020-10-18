function [PosteriorStats, PosteriorTime, Likelihood] = RunFilter2(targetData, options)


    mu0 = options.mu0;
    sigma0 = options.sigma0;
    x0 = MatToMomentVec(mu0, sigma0);
    numBins = options.NumBins;
    numMoments = options.NumMoments;
    
    tPathY = targetData.Time;
    yPath = targetData.Measurements;
    
    currT = 0;
    PosteriorStats = zeros(numMoments, (numBins - 1) * length(tPathY) - numBins);
    PosteriorTime = zeros(1, (numBins - 1) * length(tPathY) - numBins);
   
    ODEOpts.RelTol = 0.01;

    
    for k=2:length(tPathY)
        nextT = tPathY(k);

        time = linspace(0, nextT - currT, numBins);

        options.Jump = 0;
        funHandle = @(t, y) options.FilterODE(t, y, options);
        
        [~, xPrior] = dopri5Mex(funHandle, time, x0, ODEOpts);
        
        if (numBins == 2)
           xFiltered = [xPrior(1, :); xPrior(end, :)];
        else
           xFiltered = xPrior; 
        end


        options.Jump = 1;
        options.Measurement = yPath(k);
        [posteriorStat, LTmp] = PosteriorMoments(0, xFiltered(end, :), options);
        xFiltered(end, :) = xFiltered(end, :) + real(posteriorStat);
        L(k-1) = LTmp;

        x0 = xFiltered(end, :)';
        idx = (k-2)*(numBins-1)+1:(k-1)*(numBins-1);
        PosteriorStats(:, idx) = xFiltered(1:end-1, :)';
        PosteriorTime(idx) = currT + time(1:end-1);

        currT = nextT;

        if (options.StatusOutput)
            fprintf('Processed measurement %d...\n', k-2);
        end
    end

    Likelihood = sum(L);


end