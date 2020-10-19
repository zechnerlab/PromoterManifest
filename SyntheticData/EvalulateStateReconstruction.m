clear;
close all;

warning off;

addpath('../Common');
addpath('../Common/Statistics/');
addpath('../Common/StochChemKin/');
addpath('../Common/Models/');
addpath('../Common/ODEs/dopri');
addpath('../Common/ODEs/');
addpath('../Common/');
addpath('../Common/SpecialFunctions/');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');

sVec = [0, 0.1];

perturbationType = 1;

for u=1:length(sVec)
    
    perturbationStrength = sVec(u);
    
    fileName = sprintf('SyntheticDataReconstruction_%f_%d.mat', perturbationStrength, perturbationType);
    res = load(fileName);
    
    for k=1:length(res.Promoters)
        
        prom = res.Promoters{k};
        
        numCells = size(prom.MeanTauA, 1);
        
        for l=1:numCells
            
            X(l, 1) = prom.YFP{l}.totalTimeActive/60;
            X(l, 2) = prom.YFP{l}.timeToActive/60;
            X(l, 3) = prom.YFP{l}.trOutput;
            X(l, 4) = ~isnan(prom.YFP{l}.timeToActive);
            
            
        end
        
        Y(:, 1) = prom.MeanTauA/60;
        Y(:, 2) = prom.MeanTauS/60;
        Y(:, 3) = prom.MeanTrOutput;
        Y(:, 4) = prom.PActivated > 0.99;
        figure(k);
        
        responseClassAcc = sum(X(find(prom.Valid), 4) == Y(find(prom.Valid), 4)) / numCells;
        
        
        validIdx = and(~isnan(Y(:, 2)), ~isnan(X(:, 2)));
        validIdx = and(validIdx, prom.Valid);
        sigma = mean((X(validIdx, 1) - Y(validIdx, 1)).^2);

        m = fitlm(X(validIdx, 1), Y(validIdx, 1));
        coeffTauA = m.Coefficients.Estimate;
        r2 = m.Rsquared.Ordinary;

        
        subplot(length(sVec),3, (u-1)*3 + 1);
        p = plot(X(validIdx, 1), Y(validIdx, 1), 'o'); hold on;
        set(p, 'MarkerFaceColor', get(p, 'Color'));
        
        xLimits = xlim;
        xGrid = linspace(xLimits(1), xLimits(2), 20);
        yGrid = coeffTauA(1) + coeffTauA(2)*xGrid;
        yGrid2 = 0 + 1*xGrid;
        plot(xGrid, yGrid2, 'r--');
        
        xlabel('True (min)');
        ylabel('Predicted (min)');
        titleStr = sprintf('R^2=%f, k=%f', r2, coeffTauA(2));
        title(titleStr);
        
        fprintf('TimeActive, promoter %d, p=%f: d=%f,k=%f, sigma=%f, R2=%f\n', ...
            k, perturbationStrength, coeffTauA(1), coeffTauA(2), sigma, r2);
        
        validIdx = and(~isnan(Y(:, 2)), ~isnan(X(:, 2)));
        validIdx = and(validIdx, prom.Valid);
        sigma = mean((X(validIdx, 2) - Y(validIdx, 2)).^2);

        m = fitlm(X(validIdx, 1), Y(validIdx, 1));
        coeffTauS = m.Coefficients.Estimate;
        r2 = m.Rsquared.Ordinary;
        
        subplot(length(sVec),3, (u-1)*3 + 2);
        p = plot(X(validIdx, 2), Y(validIdx, 2), 'o'); hold on;
        set(p, 'MarkerFaceColor', get(p, 'Color'));
        
        xLimits = xlim;
        xGrid = linspace(xLimits(1), xLimits(2), 20);
        yGrid = coeffTauS(1) + coeffTauS(2)*xGrid;
        yGrid2 = 0 + 1*xGrid;
        plot(xGrid, yGrid2, 'r--');
        
        xlabel('True (min)');
        ylabel('Predicted (min)');
        titleStr = sprintf('R^2=%f, k=%f', r2, coeffTauS(2));
        title(titleStr);
        
        fprintf('TimeToActive, promoter %d, p=%f: d=%f,k=%f, sigma=%f, R2=%f\n', k, perturbationStrength, coeffTauA(1), coeffTauA(2), sigma, r2);        
        
        validIdx = and(~isnan(Y(:, 3)), ~isnan(X(:, 3)));
        validIdx = and(validIdx, prom.Valid);
        sigma = sqrt(mean((X(validIdx, 3) - Y(validIdx, 3)).^2));
        
        m = fitlm(X(validIdx, 1), Y(validIdx, 1));
        coeffTrO = m.Coefficients.Estimate;
        r2 = m.Rsquared.Ordinary;

        subplot(length(sVec),3, (u-1)*3 + 3);
        p = plot(X(validIdx, 3), Y(validIdx, 3), 'o'); hold on;
        set(p, 'MarkerFaceColor', get(p, 'Color'));
        
        xLimits = xlim;
        xGrid = linspace(xLimits(1), xLimits(2), 20);
        yGrid = coeffTrO(1) + coeffTrO(2)*xGrid;
        yGrid2 = 0 + 1*xGrid;
        plot(xGrid, yGrid2, 'r--');
        
        xlabel('True (transcripts)');
        ylabel('Predicted (transcripts)');
        titleStr = sprintf('R^2=%f, k=%f', r2, coeffTrO(2));
        title(titleStr);
        
        fprintf('TrOutput, promoter %d, p=%f: d=%f,k=%f, sigma=%f, R2=%f\n', k, perturbationStrength, coeffTrO(1), coeffTrO(2), sigma, r2);
        fprintf('Responder classification rate=%2.2f%%\n\n', responseClassAcc*100);
    end
end
