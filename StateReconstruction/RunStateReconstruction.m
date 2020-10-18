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
conditionsAll = GenerateConditions();

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;

%Load the promStats.mat structure, which stores for each promoter the
%maximum transcription rate across all conditions. This is needed to
%identify the promoter states which are considered "active". In particular,
%a promoter state is considered "active" if it has a transcription rate of
%at least 20% of the maximum transcription rate of that promoter. While
%this distinction is somewhat arbitrary, it provides a simple distinction
%between cells that have significant gene expression and those that don't.
res = load('../ModelCalibration/results/promStats.mat');
promStats = res.promStats;

%if set to 1, the parameter posterior samples from the preceeding MCMC run
%are used ('CalibrateModels.m');
useParameterPrior = 0;


for u=1:length(conditions)
    
    
    modelCondition = conditions(u);
    
    
    for k=1:length(promoters)
        
        
        promName = promoters{k};
        conditionName = conditions{u}.Name;
        totName = sprintf('%s_%s', promName, conditionName);
        
        %Load data and calibrated models.
        d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
        YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
        results{u} = load(['../ModelCalibration/results/ModelInference_' promName '_' modelCondition{1}.Name '.mat']);
        maxTranscription = promStats{k}.maxZ;
        
        %Use only cells which were *not* used for calibration
        cellIdx = setdiff(1:length(YFP.cells), results{u}.cellIdx);
        cells = YFP.cells(cellIdx);
        
        
        %Set the InputParams structure for each cell. Could in principle be
        %different between cells, but here it is considered to be the same.
        for j=1:length(cells)
            cells{j}.InputParams = d.TFInputParams;
            figure(99);
            plot(cells{j}.MeasurementTime/60, cells{j}.Measurement, 'o-'); hold on;
            title(promName);
        end
        drawnow;
        hold off;
        
        
        %Set model properties
        model = results{u}.modelOpt;
        config = results{u}.config;
        
        burnIn = 5000;
        posteriorSamples = results{u}.chain(:, burnIn:end);
        
        model.InputODEInfos.odeHandle = odeInfos.odeHandle;
        model.InputODEInfos.infos = odeInfos.infos;
        model.InputParams = d.TFInputParams;
        model.numBins = 3;
        model.MeasurementSigma = 0.15;

        %Filter length specifies how many time points are handeled within
        %one recursion of the Sequential Monte Carlo (SMC) algorithm. Best
        %performances seem to be achieved when only a single time point is
        %used per recursion (i.e., filterLength=2). 
        filterLength = 2;
        
        %Number of particles used for the SMC algorithm. Since mRNA and
        %protein dynamics have been integrated out, a relatively small
        %number of particles is used here.
        M = 400;
        
        %Specify sampling grid (i.e., temporal resolution of the
        %reconstructions.
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
        
        %% Run state reconstruction for each cell.
        parfor j=1:numCells
            
            %Stop time.
            tic
            cellStruct = cells{j};
            
            try
                %First initialize the particule distribution.
                if (useParameterPrior == 1)
                    particles = InitializeParticleDistribution(model, config, M, posteriorSamples);
                else
                    particles = InitializeParticleDistribution(model, config, M);
                end
                
                %Then recursively reconstruct the current cell using the hybrid SMC algorithm.
                particles = ReconstructCell(model, config, particles, cellStruct, filterLength);
                
                
                %This section post-processes the reconstructions and
                %extracts features which are used for subsequent analyses.
                %This is done on-the-fly to reduce the very large amounts
                %of generated data.
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

                    tPassed = toc;
                    
                    %Output every 10th cell to give an idea of how the
                    %algorithm progresses.
                    if (mod(j, 10)==0)
                        fprintf('Finished cell %d (tauS=%f, tauA=%f, pA=%f, n12=%f, n23=%f, numPart=%d, processing time=%f)\n', j, meanTauS/60, meanTauA/60, pActivated, meanRSwitches(1), meanRSwitches(4), numUniqueParticles, tPassed);
                    end
                else %mark cells as invalid if the SMC algorithm could not complete for numerical reasons
                    Valid(j) = 0;
                end
                
            catch err %if an error occurs, mark this cell as 'invalid'.
                Valid(j) = 0;
                particles = [];
                fprintf(err.message);
            end
            
        end
        
        
        
        %% Store all the results to the 'Promoter' structure.
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
        Promoters{k}.Valid = Valid;
        Promoters{k}.YFP = cells;
        Promoters{k}.cellIdx = cellIdx;
        Promoters{k}.NumUniqueParticles = NumUniqueParticles;
        Promoters{k}.TrReconstruction = TrReconstruction;
         
        fprintf('Processed condition %s for promoter %s (NumCells=%d, valid=%d).\n', conditionName, promName, numCells, sum(Valid));
        
    end
    
    save(['results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
end
