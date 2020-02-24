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

conditions = GenerateConditions();
for l=1:20
    condition = conditions{l};
    
    TFData = load(['../Data/Msn2/ALD3/ALD3_' condition.Name '_MSN2.mat']);
    
    
    TFInputParams = TFData.TFInputParams;
    TFInputParams.inputLevels = TFData.TFInputParams.inputLevels;
    
    config.InputParams = TFInputParams;
    
    rates{1}.Fun = @(t, u, params) params.k * u;
    rates{1}.Idx = [1, 2];
    params.k = 0.01;
    rates{1}.Params = params;
    
    
    rates{2}.Fun = @(t, u, params) params.k * (1 - u^params.n / (params.V0^params.n + u^params.n));
    rates{2}.Idx = [2, 3];
    params.k = 0.0005;
    params.n = 6;
    params.V0 = 0.12;
    rates{2}.Params = params;
    
    rates{3}.Fun = @(t, u, params) params.k;
    rates{3}.Idx = [3, 1];
    params.k = 0.003;
    rates{3}.Params = params;
    
    
    
    rates{4}.Fun = @(t, u, params) params.k *  u^params.n / (params.V0^params.n + u^params.n);
    rates{4}.Idx = [2, 4];
    params.k = 0.01;
    params.n = 6;
    params.V0 = 1.2;
    rates{4}.Params = params;
    
    rates{5}.Fun = @(t, u, params) params.k;
    rates{5}.Idx = [4, 1];
    params.k = 0.002;
    rates{5}.Params = params;
    
    
    config.rates = rates;
    
    P0 = [1; 0; 0; 0];
    
    [~, P] = ode15s(@PromoterODE, measureGrid, P0, {}, config);
    
    Z = [0, 0, 0.1, 0.1];
    
    meanTr = Z*P';
    varTr = Z.^2*P' - meanTr.^2;
    
    
    transcrOutput(l) = sum(meanTr(1:end-1).*dT);
    
    plot(measureGrid, Z*P'); hold on;
    
end



