clear all;
close all;

addpath('../../../../Common/Statistics/');
addpath('../../../../Common/StochChemKin/');
addpath('../../../../Common/Models/');
addpath('../../../../Common/ODEs/dopri');
addpath('../../../../Common/');


T = 220*60;
grid = linspace(0, T, 1000);
M = 25;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);

[Pre, Post, c, X0] = CreateThreeStageExpressionSystemLeaky();
c(1) = 0.01;
c(2) = 0.01;
c(3) = 0;
c(4) = 0.1;
c(5) = 0.001;
%c(end) = 0.00001;
X0(end) = 0;


numCells = 40;

for k=1:numCells
    [x, t] = SimulateSSA(X0, c, Pre, Post, 100000, 0, T);
    Xs = SampleCTMPPathGrid(x, t, measureGrid);
    
    measureIdx = 4;
    sigmaM = 15;
    
    protein = Xs(4, :);
    reporter = protein + sigmaM*randn(1, M);
    plot(measureGrid, reporter, 'o', t, x(end, :)); hold on;
    
    cells{k}.Measurement = reporter;
    cells{k}.MeasurementTime = measureGrid;
    cells{k}.Simulated = 1;
    cells{k}.t = t;
    cells{k}.x = x;
    cells{k}.MeasurementSigma = sigmaM;
    cells{k}.c = c;
    cells{k}.X0 = X0;
    cells{k}.Pre = Pre;
    cells{k}.Post = Post;
end


save Cells_2StatePromoter.mat cells;


