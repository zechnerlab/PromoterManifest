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
    
    fileName = sprintf('FakeDataReconstruction_%f_%d.mat', perturbationStrength, perturbationType);
    res = load(fileName);
    
    for k=1:length(res.Promoters)
        
        prom = res.Promoters{k};
        
        numCells = size(prom.MeanTauA, 1);
        
        for l=1:numCells
            
            X(l, 1) = prom.YFP{l}.totalTimeActive/60;
            X(l, 2) = prom.YFP{l}.timeToActive/60;
            X(l, 3) = prom.YFP{l}.trOutput;
            
            
        end
        
        Y(:, 1) = prom.MeanTauA/60;
        Y(:, 2) = prom.MeanTauS/60;
        Y(:, 3) = prom.MeanTrOutput;
        
        figure(k);

        
        validIdx = and(~isnan(Y(:, 1)), ~isnan(X(:, 1)));
        validIdx = and(validIdx, prom.Valid);
        sigma = mean((X(validIdx, 1) - Y(validIdx, 1)).^2);
        [coeffTauA, infoTauA] = polyfit(X(validIdx, 1), Y(validIdx, 1), 1);
        yfit = polyval(coeffTauA, X(validIdx, 1));
        resSq = sum((yfit - Y(validIdx, 1)).^2);
        resTotal = (length(Y(validIdx, 1))-1)*var(Y(validIdx, 1));
        r2 = 1 - resSq/resTotal;
        
        subplot(length(sVec),3, (u-1)*3 + 1);
        p = plot(X(validIdx, 1), Y(validIdx, 1), 'o'); hold on;
        set(p, 'MarkerFaceColor', get(p, 'Color'));
        
        xLimits = xlim;
        xGrid = linspace(xLimits(1), xLimits(2), 20);
        yGrid = coeffTauA(2) + coeffTauA(1)*xGrid;
        yGrid2 = 0 + 1*xGrid;
        plot(xGrid, yGrid2, 'r--');
        
        xlabel('True (min)');
        ylabel('Predicted (min)');
        titleStr = sprintf('R^2=%f, k=%f', r2, coeffTauA(1));
        title(titleStr);
        
        fprintf('TimeActive, promoter %d, p=%f: d=%f,k=%f, sigma=%f, R2=%f\n', k, perturbationStrength, coeffTauA(2), coeffTauA(1), sigma, r2);
        
        validIdx = and(~isnan(Y(:, 2)), ~isnan(X(:, 2)));
        validIdx = and(validIdx, prom.Valid);
        sigma = mean((X(validIdx, 2) - Y(validIdx, 2)).^2);
        [coeffTauS, infoTauS] = polyfit(X(validIdx, 2), Y(validIdx, 2), 1);
        yfit = polyval(coeffTauS, X(validIdx, 2));
        resSq = sum((yfit - Y(validIdx, 2)).^2);
        resTotal = (length(Y(validIdx, 2))-1)*var(Y(validIdx, 2));
        r2 = 1 - resSq/resTotal;
        
        subplot(length(sVec),3, (u-1)*3 + 2);
        p = plot(X(validIdx, 2), Y(validIdx, 2), 'o'); hold on;
        set(p, 'MarkerFaceColor', get(p, 'Color'));
        
        xLimits = xlim;
        xGrid = linspace(xLimits(1), xLimits(2), 20);
        yGrid = coeffTauS(2) + coeffTauS(1)*xGrid;
        yGrid2 = 0 + 1*xGrid;
        plot(xGrid, yGrid2, 'r--');
        
        xlabel('True (min)');
        ylabel('Predicted (min)');
        titleStr = sprintf('R^2=%f, k=%f', r2, coeffTauS(1));
        title(titleStr);
        
        fprintf('TimeToActive, promoter %d, p=%f: d=%f,k=%f, sigma=%f, R2=%f\n', k, perturbationStrength, coeffTauA(2), coeffTauA(1), sigma, r2);
        
        
        
        validIdx = and(~isnan(Y(:, 3)), ~isnan(X(:, 3)));
        validIdx = and(validIdx, prom.Valid);
        sigma = sqrt(mean((X(validIdx, 3) - Y(validIdx, 3)).^2));
        [coeffTrO, infoTrO] = polyfit(X(validIdx, 3), Y(validIdx, 3), 1);
        yfit = polyval(coeffTrO, X(validIdx, 3));
        resSq = sum((yfit - Y(validIdx, 3)).^2);
        resTotal = (length(Y(validIdx, 3))-1)*var(Y(validIdx, 3));
        r2 = 1 - resSq/resTotal;
        
        subplot(length(sVec),3, (u-1)*3 + 3);
        p = plot(X(validIdx, 3), Y(validIdx, 3), 'o'); hold on;
        set(p, 'MarkerFaceColor', get(p, 'Color'));
        
        xLimits = xlim;
        xGrid = linspace(xLimits(1), xLimits(2), 20);
        yGrid = coeffTrO(2) + coeffTrO(1)*xGrid;
        yGrid2 = 0 + 1*xGrid;
        plot(xGrid, yGrid2, 'r--');
        
        xlabel('True (transcripts)');
        ylabel('Predicted (transcripts)');
        titleStr = sprintf('R^2=%f, k=%f', r2, coeffTrO(1));
        title(titleStr);
        
        fprintf('TrOutput, promoter %d, p=%f: d=%f,k=%f, sigma=%f, R2=%f\n', k, perturbationStrength, coeffTrO(2), coeffTrO(1), sigma, r2);
    end
end
