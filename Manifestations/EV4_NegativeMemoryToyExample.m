clear;
close all;

addpath('../Common/');

%% Load Msn2 input data for two different pulse intervals
d{1} = load('../Data/MSN2/ALD3/ALD3_FM_4_5min_690nM_MSN2.mat');
d{2} = load('../Data/MSN2/ALD3/ALD3_FM4_20minINT_690nM_MSN2.mat');


%Define time frame
T = 150*60;
grid = linspace(0, T, 1000);
M = 200;
measureGrid = linspace(0, T, M);
measureGrid = measureGrid(1:end);


%% Initialize toy system with negative memory
[Pre, Post, c, X0] = CreateNegativeMemorySystem();
promoterSpeed = 0.1;
c(1) = 0.02;
c(2) = 0.06;
c(3) = 0.03*promoterSpeed;
c(4) = 0.2*promoterSpeed;

c(5) = 0.0006;
c(6) = 0.001;
c(7) = 0.9;
c(8) = 0.00007;

c(end) = 0.01;

ZVal = 0.01;
Z = [0, 0, ZVal, 0];
maxZ = max(Z);


%% Perform stochastic simulations of the toy promoter model for the two different Msn2 inputs
numCells = 2000;

for i=1:length(d)
    TFInputParams = d{i}.TFInputParams;
    TFInputParams.inputLevels = c(1)*d{i}.TFInputParams.inputLevels;
    
    for k=1:numCells
        
        cTmp = c;
        
        [x, t] = SimulateSSA_TV(X0, c, Pre, Post, 2000000, T, 1, TFInputParams, [0, TFInputParams.inputTimes, T], 0);
        
        Xs = SampleCTMPPathGrid_mex(x, t, measureGrid);

        Activator = Xs(5, :);
        Inhibitor = Xs(6, :);
        
        activeIdx = x(3, 1:end-1)==1;
        dT = diff(t);
        trOutput = ZVal*sum(activeIdx.*dT)*60;
        
        cells{k}.stateIdx = [1, 2, 3, 4]*Xs(1:4, :);
        cells{k}.Z = Z(cells{k}.stateIdx);
        
        ZMat(k, :) = cells{k}.Z;
        TR(i, k) = trOutput;
        
        
        if (mod(k, 100)==0)
            fprintf('Finished iteration %d (%d)\n', k, i);
        end
    end
    
    figure(1);
    subplot(1,length(d),i);
    p = plot(measureGrid/60, mean(ZMat)*60);
    
    
end

figure(2);
bar([5, 20], mean(TR'));
xlabel('Pulse interval (min)');
ylabel('Avrg. transcr. output (mol)');




