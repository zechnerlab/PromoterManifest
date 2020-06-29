clear;
close all;

warning off;

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


odeInfos = load('odeInfosLN.mat');

%measurementNoiseResults = load('MeasurementNoiseHXK1_DM_50min_3uM.mat');

conditionsAll = GenerateConditions();
%conditions = conditions(21:26);


promoters = {'HXK1', 'DCS2', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'SIP18'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;


for u=1%length(conditions)
    
    if (conditions{u}.PulseParameters(1) == 1) %% single pulse
        conc = conditions{u}.Concentration;
        modelCondition =  GetConditions(conditionsAll, [50], conc);
    else %% repeated pulses
        modelCondition =  GetConditions(conditionsAll, [50], 690);
    end
    
    for k=1:length(promoters)
        
        
        promName = promoters{k};
        conditionName = conditions{u}.Name;
        totName = sprintf('%s_%s', promName, conditionName);
        
        d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
        YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
        results{u} = load(['results/ModelInference_' promName '_' modelCondition{1}.Name '.mat']);
        
        
        if (strcmp(modelCondition{1}.Name, conditionName))
            cellIdx = setdiff(1:length(YFP.cells), results{u}.cellIdx);
        else
            cellIdx = 1:length(YFP.cells);
        end
        
        cellIdx = cellIdx(1:min(200, end));
        
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
        
        model.InputODEInfos.odeHandle = odeInfos.odeHandle;
        model.InputODEInfos.infos = odeInfos.infos;
        model.InputParams = d.TFInputParams;
        model.numBins = 3;
        model.MeasurementSigma = 0.15;%measurementNoiseResults.modelOpt.MeasurementSigma;
        fprintf('Model sigma %f\n', model.MeasurementSigma);
        
        
        %% run state reconstruction
        
        filterLength = 8;
        M = 300;
        
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
        MeanNumSwitches = zeros(numCells, 1);
        Particles = cell(numCells, 1);
        P2 = zeros(numCells, length(tGrid));
        P3 = zeros(numCells, length(tGrid));
        Mean12Switches = zeros(numCells, 1);
        Mean23Switches = zeros(numCells, 1);
        
        config.ODEType = 2;
        config.ODEStepSize = 5;
        
        for j=1:numCells
            
            tic
            cellStruct = cells{j};
            
            particles = InitializeParticleDistribution(model, M);
            particles = ReconstructCell(model, config, particles, cellStruct, filterLength);
            
            [Mu, Sigma, pActivated, P2vec, P3vec, meanTauS, meanTauA, meanNumSwitches, meanTr, meanRSwitches] = ProcessCellReconstruction(model, particles, tGrid);
            
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
            MeanTr(j, :) = meanTr;
            MeanRSwitches(j, :) = meanRSwitches;
            
            %Particles{j} = particles;
            
            tPassed = toc;
            
            if (mod(j, 10)==0)
                fprintf('Finished cell %d (tauS=%f, tauA=%f, pA=%f, n12=%f, n23=%f, processing time=%f)\n', j, meanTauS/60, meanTauA/60, pActivated, meanRSwitches(1), meanRSwitches(4), tPassed);
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
        Promoters{k}.varZ = meanZ;
        Promoters{k}.Model = model;
        Promoters{k}.tGrid = tGrid;
        Promoters{k}.P2 = P2;
        Promoters{k}.P3 = P3;
        Promoters{k}.MeanRSwitches = MeanRSwitches;
        %Promoters{k}.Particles = Particles;
        %Promoters{k}.Valid = Valid;
        Promoters{k}.YFP = cells;
        Promoters{k}.cellIdx = cellIdx;
        
        
        fprintf('Processed condition %s for promoter %s.\n', conditionName, promName);
        
    end
    
    legend(legendStr);
    
    save(['results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
end
