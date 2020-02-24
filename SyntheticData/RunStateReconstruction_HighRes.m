clear;
close all;

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


warning off;

odeInfos = load('odeInfosLN.mat');

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

speedVec = [1, 10];

numConsideredCells = 3;



for u=1:length(speedVec)
    
    fileName = sprintf('../Data/Promoter_s%d_fakeData.mat', speedVec(u));
    load(fileName);
    
    for j=1:numConsideredCells
        usedCells{j} = cells{j};
        usedCells{j}.InputParams = model.InputParams;
        figure(99);
        plot(usedCells{j}.MeasurementTime/60, usedCells{j}.Measurement, 'o-'); hold on;
    end
    drawnow;
    hold off;
    
    
    model.InputODEInfos.odeHandle = odeInfos.odeHandle;
    model.InputODEInfos.infos = odeInfos.infos;
    model.NumBins = 3;
    
    %% run state reconstruction
    
    filterLength = 15;
    M = 1000;
    
    numCells = length(usedCells);
    tGrid = linspace(0, max(usedCells{1}.MeasurementTime), 100);
    
    meanRna = zeros(numCells, length(tGrid));
    meanProtein = zeros(numCells, length(tGrid));
    meanZ = zeros(numCells, length(tGrid));
    varRna = zeros(numCells, length(tGrid));
    varProtein = zeros(numCells, length(tGrid));
    varZ = zeros(numCells, length(tGrid));
    
    MeanTauS = zeros(numCells, 1);
    MeanTauA = zeros(numCells, 1);
    PActivated = zeros(numCells, 1);
    MeanTr = zeros(numCells, length(tGrid));
    Q95Tr = zeros(numCells, length(tGrid));
    Q5Tr = zeros(numCells, length(tGrid));
    VarTr = zeros(numCells, length(tGrid));
    MeanNumSwitches = zeros(numCells, 1);
    Particles = cell(numCells, 1);
    P2 = zeros(numCells, length(tGrid));
    P3 = zeros(numCells, length(tGrid));
    Mean12Switches = zeros(numCells, 1);
    Mean23Switches = zeros(numCells, 1);
    MeanRSwitches = zeros(numCells, 6);
    MeanTauState = zeros(numCells, 3);
    MeanTrOutput = zeros(numCells, 1);
    Valid = zeros(numCells, 1);
    
    config.ODEType = 2;
    config.ODEStepSize = 20;
    
    
    
    for j=1:numCells
        
        modelTmp = model;
        
        tic
        cellStruct = usedCells{j};
        cellIdx(j) = j;
        
        particles = InitializeParticleDistribution(modelTmp, config, M);
        particles = ReconstructCell(modelTmp, config, particles, cellStruct, filterLength);
        
        
        if (~isempty(particles))
            [Mu, Sigma, pActivated, P2vec, P3vec, meanTauS, meanTauA, meanNumSwitches, meanTr, varTr, q5Tr, q95Tr, meanTrOutput, meanRSwitches, meanTauState, ZMat] = ...
                ProcessCellReconstruction(modelTmp, particles, tGrid, max(modelTmp.Z));
            
            meanRna(j, :) = Mu(1, :);
            meanProtein(j, :) = Mu(2, :);
            meanZ(j, :) = Mu(3, :);
            
            varRna(j, :) = Sigma(1, :);
            varProtein(j, :) = Sigma(2, :);
            varZ(j, :) = Sigma(3, :);
            
            P2(j, :) = P2vec;
            P3(j, :) = P3vec;
            
            MeanTauS(j) = meanTauS;
            MeanTauA(j) = meanTauA;
            PActivated(j) = pActivated;
            
            MeanNumSwitches(j) = meanNumSwitches;
            MeanTauState(j, :) = meanTauState;
            MeanTr(j, :) = meanTr;
            VarTr(j, :) = varTr;
            Q95Tr(j, :) = q95Tr;
            Q5Tr(j, :) = q5Tr;
            
            MeanRSwitches(j, :) = meanRSwitches;
            Valid(j) = 1;
            MeanTrOutput(j) = meanTrOutput;
            
            
            tPassed = toc;
            
            if (mod(j, 1)==0)
                fprintf('Finished cell %d (tauS=%f, tauA=%f, pA=%f, n12=%f, n23=%f, processing time=%f)\n', j, meanTauS/60, meanTauA/60, pActivated, meanRSwitches(1), meanRSwitches(4), tPassed);
            end
            
        else
            Valid(j) = 0;
        end
        
        
        %plot(tGrid, meanProtein(j, :), tGrid, meanProtein(j, :) - sqrt(varProtein(j, :)), tGrid, meanProtein(j, :) + sqrt(varProtein(j, :)), cells{j}.MeasurementTime, cells{j}.Measurement, 'o'); drawnow;
        
    end
    
    
    
    %% storing
    
    Promoters{u}.meanRna = meanRna;
    Promoters{u}.meanProtein = meanProtein;
    Promoters{u}.meanZ = meanZ;
    Promoters{u}.varRna = varRna;
    Promoters{u}.varProtein = varProtein;
    Promoters{u}.PActivated = PActivated;
    Promoters{u}.MeanTauS = MeanTauS;
    Promoters{u}.MeanTauA = MeanTauA;
    Promoters{u}.MeanNumSwitches = MeanNumSwitches;
    Promoters{u}.MeanTr = MeanTr;
    Promoters{u}.VarTr = VarTr;
    Promoters{u}.Q95Tr = Q95Tr;
    Promoters{u}.Q5Tr = Q5Tr;
    Promoters{u}.MeanTrOutput = MeanTrOutput;
    Promoters{u}.MeanTauState = MeanTauState;
    Promoters{u}.varZ = meanZ;
    Promoters{u}.Model = model;
    Promoters{u}.tGrid = tGrid;
    Promoters{u}.P2 = P2;
    Promoters{u}.P3 = P3;
    Promoters{u}.MeanRSwitches = MeanRSwitches;
    Promoters{u}.Valid = Valid;
    Promoters{u}.YFP = cells;
    Promoters{u}.cellIdx = cellIdx;
    
end


fileName = sprintf('FakeDataReconstruction_HighRes.mat');
save(fileName);

