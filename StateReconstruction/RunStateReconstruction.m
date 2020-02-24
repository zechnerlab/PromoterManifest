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

%parpool(35);

warning off;

odeInfos = load('odeInfosLN.mat');

%measurementNoiseResults = load('MeasurementNoiseHXK1_DM_50min_3uM.mat');

conditionsAll = GenerateConditions();
%conditions = conditions(21:26);


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'SIP18'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;

res = load('../ModelCalibration/results/promStats.mat');
promStats = res.promStats;

for u=1:30%length(conditions)
    
    if (conditions{u}.PulseParameters(1) == 1) %% single pulse
        conc = conditions{u}.Concentration;
        modelCondition =  GetConditions(conditionsAll, [50], conc);
        
    else %% repeated pulses
        modelCondition =  conditions(21);%GetConditions(conditionsAll, [10], 690);
    end
    
    modelCondition = conditions(u);
    
    
    for k=1:length(promoters)
        
        
        promName = promoters{k};
        conditionName = conditions{u}.Name;
        totName = sprintf('%s_%s', promName, conditionName);
        
        d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
        YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
        results{u} = load(['../ModelCalibration/results/ModelInference_' promName '_' modelCondition{1}.Name '.mat']);
        maxTranscription = promStats{k}.maxZ;
        
        if (strcmp(modelCondition{1}.Name, conditionName))
            cellIdx = setdiff(1:length(YFP.cells), results{u}.cellIdx);
        else
            cellIdx = 1:length(YFP.cells);
        end
        
        %cellIdx = cellIdx(1:min(200, end));
        
        cells = YFP.cells(cellIdx);
        
        for j=1:length(cells)
            cells{j}.InputParams = d.TFInputParams;
            figure(99);
            plot(cells{j}.MeasurementTime/60, cells{j}.Measurement, 'o-'); hold on;
            title(promName);
        end
        drawnow;
        hold off;
        
        
        model = results{u}.modelOpt;
        burnIn = 5000;
        posteriorSamples = results{u}.chain(:, burnIn:end);
        config = results{u}.config;
        model.InputODEInfos.odeHandle = odeInfos.odeHandle;
        model.InputODEInfos.infos = odeInfos.infos;
        model.InputParams = d.TFInputParams;
        model.numBins = 3;
        model.MeasurementSigma = 0.15;

        
        %% run state reconstruction
        
        filterLength = 2;
        M = 400;
        
        tGrid = linspace(0, max(YFP.cells{1}.MeasurementTime), 100);
        
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
        TrReconstruction = cell(numCells, 1);
        NumUniqueParticles = zeros(numCells, 1);
        
        
        config.ODEType = 2;
        config.ODEStepSize = 20;
        
        config.UnknownCIdx = [];
        config.UnknownHIdx = [];
        config.UnknownMIdx = [];
        
        parfor j=1:numCells
            
            tic
            cellStruct = cells{j};
            
            try
                %particles = InitializeParticleDistribution(model, config, M, posteriorSamples);
                particles = InitializeParticleDistribution(model, config, M);
                particles = ReconstructCell(model, config, particles, cellStruct, filterLength);
                
                
                if (~isempty(particles))
                    [Mu, Sigma, pActivated, P2vec, P3vec, meanTauS, meanTauA, meanNumSwitches, meanTr, varTr, meanTrOutput, meanRSwitches, meanTauState, numUniqueParticles, MAPIdx] = ...
                        ProcessCellReconstruction(model, particles, tGrid, maxTranscription);
                    
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
                    
                    NumUniqueParticles(j) = numUniqueParticles;
                    
                    TrReconstruction{j}.t = particles{MAPIdx}.t;
                    TrReconstruction{j}.Z = particles{MAPIdx}.Z;
                    
                    
                    %                 for i=1:length(particles)
                    %                     particlesR{i}.t = particles{i}.t;
                    %                     particlesR{i}.Z = particles{i}.Z;
                    %                     particlesR{i}.stateIdx = particles{i}.stateIdx;
                    %                     particlesR{i}.r = particles{i}.r;
                    %                     particlesR{i}.Idx = particles{i}.Idx;
                    %                 end
                    %Particles{j} = particlesR;
                    
                    tPassed = toc;
                    
                    if (mod(j, 10)==0)
                        fprintf('Finished cell %d (tauS=%f, tauA=%f, pA=%f, n12=%f, n23=%f, numPart=%d, processing time=%f)\n', j, meanTauS/60, meanTauA/60, pActivated, meanRSwitches(1), meanRSwitches(4), numUniqueParticles, tPassed);
                    end
                    
                    
                else
                    Valid(j) = 0;
                end
                
            catch err
                Valid(j) = 0;
                particles = [];
                %Particles{j} = {};
                fprintf(err.message);
                
            end
            
            
            %plot(tGrid, meanProtein(j, :), tGrid, meanProtein(j, :) - sqrt(varProtein(j, :)), tGrid, meanProtein(j, :) + sqrt(varProtein(j, :)), cells{j}.MeasurementTime, cells{j}.Measurement, 'o'); drawnow;
            
        end
        
        
        
        %% storing
        
        
        Promoters{k}.meanRna = meanRna;
        Promoters{k}.meanProtein = meanProtein;
        Promoters{k}.meanZ = meanZ;
        Promoters{k}.varRna = varRna;
        Promoters{k}.varProtein = varProtein;
        Promoters{k}.PActivated = PActivated;
        Promoters{k}.MeanTauS = MeanTauS;
        Promoters{k}.MeanTauA = MeanTauA;
        Promoters{k}.MeanNumSwitches = MeanNumSwitches;
        Promoters{k}.MeanTr = MeanTr;
        Promoters{k}.VarTr = VarTr;
        Promoters{k}.MeanTrOutput = MeanTrOutput;
        Promoters{k}.MeanTauState = MeanTauState;
        Promoters{k}.varZ = meanZ;
        Promoters{k}.Model = model;
        Promoters{k}.tGrid = tGrid;
        Promoters{k}.P2 = P2;
        Promoters{k}.P3 = P3;
        Promoters{k}.MeanRSwitches = MeanRSwitches;
        %Promoters{k}.Particles = Particles;
        Promoters{k}.Valid = Valid;
        Promoters{k}.YFP = cells;
        Promoters{k}.cellIdx = cellIdx;
        Promoters{k}.NumUniqueParticles = NumUniqueParticles;
        Promoters{k}.TrReconstruction = TrReconstruction;
        
        
        fprintf('Processed condition %s for promoter %s (NumCells=%d, valid=%d).\n', conditionName, promName, numCells, sum(Valid));
        
    end
    
    legend(legendStr);
    
    save(['results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
end
