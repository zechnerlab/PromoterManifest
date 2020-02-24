clear;
close all;

addpath('../Common/');
addpath('../../../Common/Statistics');
warning off;

%d{1} = load('../Data/MSN2/ALD3/ALD3_FM_4_5min_690nM_MSN2.mat');
d{1} = load('../Data/MSN2/ALD3/ALD3_DM_50min_275nM_MSN2.mat');
d{2} = load('../Data/MSN2/ALD3/ALD3_FM4_20minINT_690nM_MSN2.mat');



T = 150*60;
grid = linspace(0, T, 1000);
M = 200;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);


promoterSpeed = 0.1;
[Pre, Post, c, X0] = CreateComplexPromoterModel();


numCells = 2000;

TFInputParams = d{1}.TFInputParams;
TFInputParams.inputLevels = c(1)*d{1}.TFInputParams.inputLevels;


cIdx = [1, 5];
cFun{1} = @(t, X, params) X(1) * GetPieceWiseConstantInput(t, TFInputParams);
cFun{2} = @(t, X, params) X(2) * params.k * X(4)^params.n / (params.V0^params.n + X(4)^params.n);
params.k = 0.01;
params.n = 50;
params.V0 = 10;
cFunParams{1} = {};%params;
cFunParams{2} = params;

[x, t] = SimulateGillespieRL(X0, c, Pre, Post, 100000, T, cIdx, cFun, cFunParams);
Xs = SampleCTMPPathGrid_mex(x, t, measureGrid);

stairs(t, x(4, :)); hold on;
stairs(t, x(3, :));


