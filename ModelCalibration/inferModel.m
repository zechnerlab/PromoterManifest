function Results = inferModel(promName, conditionName, config, K, firstModel)


        targetData = struct;
        totName = sprintf('%s_%s', promName, conditionName);
        
        d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
        YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
        
        rndIdx = randperm(length(YFP.cells));
        
        %use half of the cells for calibration (and the other half for
        %state reconstruction + validation).
        numCells = floor(0.5*length(rndIdx));
        cellIdx = rndIdx(1:numCells);
        tmpCells = YFP.cells(cellIdx);
        
        Meas0 = [];
        
        Y = zeros(numCells, length(tmpCells{1}.Measurement));
        for i=1:numCells
            
            Y(i, :) = tmpCells{i}.Measurement;
            Meas0(end+1) = tmpCells{i}.Measurement(1);
            
            tmpCells{i}.InputParams = d.TFInputParams;
        end
        
        
        data = cell(size(tmpCells{1}.Measurement));
        for i=1:length(tmpCells{i}.Measurement)
            data{i} = Y(:, i);
        end
        
        
        %Allows to exclude time-points from inference. Here, all time
        %points are used.
        validIdx = 1:length(tmpCells{1}.Measurement);
        
        %Determine moment uncertainties via bootstrapping
        [Moments, Uncertainties] = BootstrapMomentsStd(data, 2000);
        
        
        targetData.Conditions{1}.Moments = Moments;
        targetData.Conditions{1}.Uncertainties = Uncertainties;
        targetData.Conditions{1}.InputParams = d.TFInputParams;
        targetData.Conditions{1}.CellIdx = cellIdx;
        targetData.Time = tmpCells{1}.MeasurementTime(1:end);
        
    
        %Determine initial mean and variance of YFP abundance.
        muMeas0 = mean(Meas0);
        varMeas0 = var(Meas0);
        
        %Either a model is provided as a starting point, or it is created
        %using default parameters.
        if (nargin<5)
            firstModel = CreateDefaultModel(3, config, d.TFInputParams, muMeas0, varMeas0);
        end
        
        model = firstModel;
        model.muMeas0 = muMeas0;
        model.varMeas0 = varMeas0;
        
        %Determine start parameters for inference.
        CommonParams = GetModelParameters(model, config);
        startParams = CommonParams;
        
        %Perform MCMC parameter inference. Returns MAP estimates.
        [optParams, ~, chain] = RunMCMCMoments(model, config, targetData, startParams, K, 0.05, validIdx);
        
        
        %Store best-fitting parameters into a model struture.
        modelOpt = SetModelParameters(model, config, optParams);
        
        
        Results.model = model;
        Results.modelOpt = modelOpt;
        Results.chain = chain;
        Results.targetData = targetData;
        Results.cellIdx = cellIdx;
        Results.config = config;
        Results.K = K;
        Results.validIdx = validIdx;
        Results.promName = promName;
        Results.startParams = startParams;
        
        
        
        