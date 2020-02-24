clear;
close all;

addpath('../Common/');
addpath('../../../Common/Statistics');
addpath('../../../Common/Models');
addpath('results/');
addpath('../Data/MSN2');

warning off;

T = 80*60;
grid = linspace(0, T, 1000);
M = 200;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);

twoState = 1;

%% Set up model
[Pre, Post, c, X0] = CreateThreeStateExpressionSystemLeaky();

trMaxRate = 0.5;

if twoState
    %two-state promoters
    c(1) = 0.001;
    c(2) = 0.001;
    c(3) = 0;
    c(4) = 0;
    c(5) = trMaxRate;
    c(6) = 0;
    
else
    %three-state promoters
    c(1) = 0.01;
    c(2) = 0.01;
    c(3) = 0.006;
    c(4) = 0.001;
    c(5) = 0.1*trMaxRate;
    c(6) = trMaxRate;
    c([8, 9]) = 0;
    c(10) = 0.001;
    
end

Z = [0, c(5), c(6)];
gActIdx = find(Z>0.2*max(Z));


cVec = [100];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();

conditions = GetConditions(conditions, plVec, cVec);

for u=1:length(conditions)
    condition = conditions{u};
    promoters = {'ALD3'};
    TFData = load(['../Data/Msn2/ALD3/ALD3_' condition.Name '_MSN2.mat']);
    
    
    TFInputParams = TFData.TFInputParams;
    TFInputParams.inputLevels = TFData.TFInputParams.inputLevels;
    
    
    cIdx = [1];
    cFun{1} = @(t, X, c1) X(1) * c1 * GetPieceWiseConstantInput(t, TFInputParams);
    cFunParams{1} = c(1);
    
    
    for k=1:50
        [x, t, r] = SimulateGillespieRL(X0, c, Pre, Post, 100000, T, cIdx, cFun, cFunParams);
        Xs = SampleCTMPPathGrid_mex(x, t, measureGrid);
        
        rnaIdx = 4;
        transcrIdx = [5, 6];
        txEvents(k, u) = sum(sum(repmat(r, 2, 1)==transcrIdx'));
        

        gAct = sum(x(gActIdx, :)==1, 1);
        tauSIdx = find(gAct);
        
        if (sum(tauSIdx)>0)
        
        dT = t(2:end) - t(1:end-1);
        tauA = dT*gAct(1:end-1)';
        
        timeToActivate(k, u) = t(tauSIdx(1));
        timeActive(k, u) = tauA;
        responder(k, u) = 1;
        else
           timeToActivate(k, u) = nan;
           timeActive(k, u) = nan;
           responder(k, u) = 0; 
        end
        
        ZVec(k, :) = Z*Xs([1, 2, 3], :);
    end
    
    plMat(u) = condition.PulseParameters(2);
    concMat(u) = condition.Concentration;
    maxTr(u) = max(mean(ZVec));
    
    fprintf('Finished condition %d (%d)\n', u, length(conditions));
    
end

subplot(1,3,1);
loglog(mean(txEvents), var(txEvents)./mean(txEvents).^2, 'o'); hold on;

subplot(1,3,2);
plot(maxTr, nanmean(timeToActivate), 'o');

subplot(1,3,3);
plot(maxTr, nanmean(timeActive), 'o');

figure;
plot(plMat, mean(txEvents), 'o');


