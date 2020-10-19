%Runs trajectory inference using synthetic data.

clear;
close all;

addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../StateReconstruction/');
addpath('../');


warning off;

odeInfos = load('odeInfosLN.mat');

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

speedVec = [1, 10];
pVec = [0, 0.1];
perturbParams = 1;

for o=1:length(pVec)
    perturbationStrength = pVec(o);
    
    for u=1:length(speedVec)
        
        fileName = sprintf('../Data/Promoter_s%d_synData.mat', speedVec(u));
        load(fileName);
        
        
        for j=1:length(cells)
            cells{j}.InputParams = model.InputParams;
            figure(99);
            plot(cells{j}.MeasurementTime/60, cells{j}.Measurement, 'o-'); hold on;
        end  
        drawnow;
        hold off;
        
        
        model.InputODEInfos.odeHandle = odeInfos.odeHandle;
        model.InputODEInfos.infos = odeInfos.infos;
        model.NumBins = 3;
        model.MeasurementSigma = 0.05;
        %% run state reconstruction
        
        filterLength = 2;
        M = 400;
        
        tGrid = linspace(0, max(cells{1}.MeasurementTime), 100);
        
        numCells = length(cells);
        
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
        Models = cell(numCells, 1);
        
        config.ODEType = 2;
        config.ODEStepSize = 20;
        
        modelPert = model;
        modelPert = rmfield(modelPert, 'cells');
        
        if (perturbParams == 2)
            params = model.Q(:);
            paramsPerturbed = lognrnd(log(params), perturbationStrength);
            Q = reshape(paramsPerturbed, 3, 3);
            Q = Q - diag(diag(Q));
            Q = Q - diag(sum(Q));
            modelPert.Q = Q;
            
            modelPert.Z(2:3) = lognrnd(log(modelPert.Z(2:3)), perturbationStrength);
            %modelPert.c = lognrnd(log(modelPert.c), perturbationStrength);
            %modelPert.MeasurementSigma = lognrnd(log(modelPert.MeasurementSigma), perturbationStrength);
        end
        
        parfor j=1:numCells
            
            modelTmp = modelPert;
            if (perturbParams == 1)
                
                params = model.Q(:);
                paramsPerturbed = lognrnd(log(params), perturbationStrength);
                Q = reshape(paramsPerturbed, 3, 3);
                Q = Q - diag(diag(Q));
                Q = Q - diag(sum(Q));
                modelTmp.Q = Q;

                modelTmp.Z(2:3) = lognrnd(log(modelTmp.Z(2:3)), perturbationStrength);
                modelTmp.c = lognrnd(log(modelTmp.c), perturbationStrength);
                modelTmp.MeasurementSigma = lognrnd(log(modelTmp.MeasurementSigma), perturbationStrength);
            end
             
            tic
            cellStruct = cells{j};
            cellIdx(j) = j;
            
            particles = InitializeParticleDistribution(modelTmp, config, M);
            particles = ReconstructCell(modelTmp, config, particles, cellStruct, filterLength);
            
            
            if (~isempty(particles))
                [Mu, Sigma, pActivated, P2vec, P3vec, meanTauS, meanTauA, meanNumSwitches, meanTr, varTr, meanTrOutput, meanRSwitches, meanTauState, numUniqueParticles, MAPIdx] = ...
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
                MeanRSwitches(j, :) = meanRSwitches;
                Valid(j) = 1;
                MeanTrOutput(j) = meanTrOutput;
                Models{j} = modelTmp;

                tPassed = toc;
                
                if (mod(j, 1)==0)
                    fprintf('Finished cell %d (tauS=%f, tauA=%f, pA=%f, n12=%f, n23=%f, processing time=%f)\n', j, meanTauS/60, meanTauA/60, pActivated, meanRSwitches(1), meanRSwitches(4), tPassed);
                end
                
            else
                Valid(j) = 0;
            end

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
        Promoters{u}.MeanTrOutput = MeanTrOutput;
        Promoters{u}.MeanTauState = MeanTauState;
        Promoters{u}.varZ = meanZ;
        Promoters{u}.Model = model;
        Promoters{u}.tGrid = tGrid;
        Promoters{u}.P2 = P2;
        Promoters{u}.P3 = P3;
        Promoters{u}.MeanRSwitches = MeanRSwitches;
        %Promoters{k}.Particles = Particles;
        Promoters{u}.Valid = Valid;
        Promoters{u}.YFP = cells;
        Promoters{u}.cellIdx = cellIdx;
        Promoters{u}.Models = Models;
    end
    
    
    fileName = sprintf('SyntheticDataReconstruction_%f_%d.mat', perturbationStrength, perturbParams);
    save(fileName);
end
