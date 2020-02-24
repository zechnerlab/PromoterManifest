clear;
close all;

addpath('../../Common/Statistics/');
addpath('../../Common/StochChemKin/');
addpath('../../Common/Models/');
addpath('../../Common/ODEs/dopri');
addpath('../../Common/');
addpath('../../Common/');
addpath('../Common/');

d = load('MSN2/ALD3/ALD3_DM_40min_690nM_MSN2.mat');



T = 150*60;
grid = linspace(0, T, 1000);
M = 55;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);


promoterSpeed = 1;
[Pre, Post, c, X0] = CreateThreeStateExpressionSystemLeaky();
c(1) = 0.05;
c(2) = 0.055;
c(3) = 0.001*promoterSpeed;
c(4) = 0.004*promoterSpeed;

c(5) = 0.0035;
c(6) = 0.728;
c(7) = 0.0013;
c(8) = 0.1;
c(9) = 1.6667e-05;
%c(end) = 0.00005;
%X0(end) = 0;

Z = [0, c(5), c(6)];
maxZ = max(Z);

activeState = find(Z>0.2*maxZ);

eIdx = 8;
extrSigmaSquared = 0.02^2;
extrMean = 1;%c(eIdx);
extrCV = sqrt(extrSigmaSquared) / extrMean;
a = 1/extrCV^2;
b = a / extrMean;

numCells = 30;

TFInputParams = d.TFInputParams;
TFInputParams.inputLevels = c(1)*d.TFInputParams.inputLevels;

for k=1:numCells
    
    cTmp = c;
    cTmp(eIdx) = c(eIdx)*gamrnd(a, 1/b);
    
    [x, t] = SimulateSSA_TV(X0, c, Pre, Post, 2000000, T, 1, TFInputParams, [0, TFInputParams.inputTimes, T], 0);
    
    %[x, t] = SimulateSSA(X0, cTmp, Pre, Post, 100000, 0, T);
    Xs = SampleCTMPPathGrid(x, t, measureGrid);
    
    measureIdx = 4;
    sigmaM = 0.05;
    
    protein = Xs(5, :);
    reporter = lognrnd(log(protein+eps), sigmaM); %needed to avoid negative values
    subplot(1,2,1);
    plot(measureGrid, reporter, 'o', t, x(5, :)); hold on;
    
    subplot(1,2,2);
    plot(t, x(4, :)); hold on;
    
    cells{k}.Measurement = reporter;
    cells{k}.MeasurementTime = measureGrid;
    cells{k}.Simulated = 1;
    cells{k}.t = t;
    cells{k}.x = x;
    cells{k}.stateIdx = [1, 2, 3]*x(1:3, :);
    cells{k}.Z = Z(cells{k}.stateIdx);
    
    dT = diff(cells{k}.t);
    for j=1:3
        tauState(k, j) = sum(dT(find(cells{k}.stateIdx(1:end-1)==j)));
    end
    cells{k}.tauState = tauState;
    
    cells{k}.trOutput = tauState(k, :)*Z';
    
    
    S = [];
    
    for l=1:length(activeState)
        S(l, :) =  cells{k}.stateIdx == activeState(l);
    end
    
    responders = sum(S, 1)>0;
    numSwitches(k) = sum(and(responders(1:end-1)==0, responders(2:end)==1));
    activeIdx = find(responders);
    if (isempty(activeIdx))
        timeToActive = nan;
    else
        timeToActive = cells{k}.t(activeIdx(1));
    end
    
    deltaT = [diff(cells{k}.t), 0];
    
    cells{k}.totalTimeActive = sum(deltaT(activeIdx));
    cells{k}.numSwitches = numSwitches;
    cells{k}.timeToActive = timeToActive;
    
end

model.c = [1; c(7:end-1)];
model.InputParams = d.TFInputParams;
model.H(1) = a/b;
model.H(2) = a/b^2;
model.muMeas0 = 0;
model.varMeas0 = 0;
model.cells = cells;
model.MeasurementSigma = sigmaM;
model.X0 = X0;
model.Pre = Pre;
model.Post = Post;
model.P0 = [1;0;0];
model.Z = Z; 
model.Q = [ 0       c(2)    0
            c(1)    0       c(4)
            0       c(3)    0 ];
model.Q = model.Q - diag(sum(model.Q));
model.StatusOutput = 0;
%Data.Model = model;

paramStr = PrintModelParams(model);

fprintf('%s\n', paramStr);

fileName = sprintf('Promoter_s%d_fakeData.mat', promoterSpeed);
save(fileName);


