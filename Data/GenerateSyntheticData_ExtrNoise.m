clear all;
close all;

addpath('../../../../Common/Statistics/');
addpath('../../../../Common/StochChemKin/');
addpath('../../../../Common/Models/');
addpath('../../../../Common/ODEs/dopri');
addpath('../../../../Common/');


T = 220*60;
grid = linspace(0, T, 1000);
M = 35;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);

[Pre, Post, c, X0] = CreateThreeStageExpressionSystemLeaky();
c(1) = 0.02;
c(2) = 0.01;
c(3) = 0;
c(4) = 0.1;
c(5) = 0.001;
%c(end) = 0.00001;
X0(end) = 0;

eIdx = 6;
extrCV = 0.4;
extrMean = c(eIdx);
a = 1/extrCV^2;
b = a / extrMean;

numCells = 100;

for k=1:numCells
    
    cTmp = c;
    cTmp(eIdx) = gamrnd(a, 1/b);
    
    [x, t] = SimulateSSA(X0, cTmp, Pre, Post, 100000, 0, T);
    Xs = SampleCTMPPathGrid(x, t, measureGrid);
    
    measureIdx = 4;
    sigmaM = 15;
    
    protein = Xs(4, :);
    reporter = max(0, protein + sigmaM*randn(1, M)); %needed to avoid negative values
    plot(measureGrid, reporter, 'o', t, x(end, :)); hold on;
    
    cells{k}.Measurement = reporter;
    cells{k}.MeasurementTime = measureGrid;
    cells{k}.Simulated = 1;
    cells{k}.t = t;
    cells{k}.x = x;
    cells{k}.c = cTmp;
    cells{k}.c(eIdx) = 1; %set to one because drawn from distribution


end

Data.H(1) = a/b;
Data.H(2) = a/b^2;
Data.cells = cells;
Data.MeasurementSigma = sigmaM;
Data.X0 = X0;
Data.Pre = Pre;
Data.Post = Post;

save Cells_2StatePromoter.mat Data;


