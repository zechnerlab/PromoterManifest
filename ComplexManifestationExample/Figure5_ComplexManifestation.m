clear;
close all;

addpath('../Common/');
addpath('../../../Common/Statistics');
addpath('../../../Common/Models');
addpath('results/');
addpath('../Data/MSN2');

warning off;

T = 120*60;
grid = linspace(0, T, 1000);
M = 200;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);
dT = measureGrid(2:end) - measureGrid(1:end-1);

concVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];

config.numPromoterStates = 4;
config.Z = [0, 0.3, 0, 0.8];

conditions = GenerateConditions();
numRows = 5;
numCols = 6;

for l=1:30
    condition = conditions{l};
    
    TFData = load(['../Data/Msn2/ALD3/ALD3_' condition.Name '_MSN2.mat']);
    
    
    TFInputParams = TFData.TFInputParams;
    TFInputParams.inputLevels = TFData.TFInputParams.inputLevels;
    
    config.InputParams = TFInputParams;
    
    rates{1}.Fun = @(t, u, params) params.k * u;
    rates{1}.Idx = [1, 2];
    params.k = 0.01;
    rates{1}.Params = params;
    
    rates{5}.Fun = @(t, u, params) params.k;
    rates{5}.Idx = [2, 1];
    params.k = 0.01;
    rates{5}.Params = params;
    
    rates{2}.Fun = @(t, u, params) params.k * (1 - u^params.n / (params.V0^params.n + u^params.n));
    rates{2}.Idx = [2, 3];
    params.k = 0.01;
    params.n = 6;
    params.V0 = 0.5;
    rates{2}.Params = params;
    
    rates{3}.Fun = @(t, u, params) params.k * (u^params.n / (params.V0^params.n + u^params.n));
    rates{3}.Idx = [3, 2];
    params.k = 0.01;
    params.n = 2;
    params.V0 = 0.001;
    rates{3}.Params = params;
    
    rates{4}.Fun = @(t, u, params) params.k *  u^params.n / (params.V0^params.n + u^params.n);
    rates{4}.Idx = [2, 4];
    params.k = 0.1;
    params.n = 3;
    params.V0 = 1.2;
    rates{4}.Params = params;
    
    rates{5}.Fun = @(t, u, params) params.k;
    rates{5}.Idx = [4, 1];
    params.k = 0.01;
    rates{5}.Params = params;
    
    
    config.rates = rates;
    
    P0 = [1; 0; 0; 0];
    M0 = [0; 0; 0; 0; 0; 0; 0; 0;];
    
    [~, X] = ode15s(@PromoterODE, measureGrid, [P0(:); M0(:)], {}, config);
    P = X(:, 1:config.numPromoterStates);
    M = X(:, config.numPromoterStates+1:end);
    
    condM = M(:, 1:config.numPromoterStates);
    condS = M(:, config.numPromoterStates+1:end);
    
    
    meanTr = config.Z*P';
    varTr = config.Z.^2*P' - meanTr.^2;
    
    
    meanR = sum(condM, 2);
    varR = sum(condS, 2) - meanR.^2;
    
    
    transcrOutput(l) = sum(meanTr(1:end-1).*dT);
    
    meanTrI = config.Z.*P;
    transcrOutputI(l, :) = sum(meanTrI(1:end-1, :).*dT');
    
    if (condition.PulseParameters(1) == 1)
       concIdx = find(concVec == condition.Concentration);
       plIdx = find(plVec == condition.PulseParameters(2));
       PulseTr(concIdx, plIdx) = meanR(end); 
       PulseTrVar(concIdx, plIdx) = varR(end);
    end
    
    QNet = CalculateFlows(P, measureGrid, config);
    QNet = QNet - diag(diag(QNet));
    
%     for u=1:config.numPromoterStates
%        QNet(u, u) = nan; 
%     end
    
    %plot(measureGrid, config.Z*P'); hold on;
    
    %plot(measureGrid, meanR, measureGrid, meanR - sqrt(varR), measureGrid, meanR + sqrt(varR)); hold on;
    subplot(numRows, numCols, l);
    heatmap(QNet, 'CellLabelColor', 'none', 'GridVisible', 'off', 'ColorbarVisible', 'off');
    title(condition.Name);
    
    %image(QNet);
end

%loglog(PulseTr, PulseTrVar./PulseTr.^2, 'o');


