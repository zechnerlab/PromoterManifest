clear;
close all;

addpath('../Data/Msn2');

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};


legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();




rows = 5;
cols = 6;

stdFactor = 1;

for u=1:30
    
    
    for k=1:length(promoters)
        promName = promoters{k};
        
        
        totName = sprintf('%s_%s', promName, conditions{u}.Name);
        YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
        res = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
        
        concIdx = find(conditions{u}.Concentration == cVec);
        plIdx = floor(conditions{u}.PulseParameters(2)/10);
        
        prom = res.Promoters{k};
        model = prom.Model;
        time = prom.YFP{1}.MeasurementTime;
        
        
        M0 = model.InputODEInfos.infos.DefaultInitialConditions;
        
        model.InitialMeans = 0;
        model.InitialVars = 0;
        means = [model.InitialMeans; model.muMeas0; model.H(1)];
        vars = [model.InitialVars; model.varMeas0; model.H(2)];
        
        M0 = SetInitialInputMoments(model.InputODEInfos.infos, M0, means, vars);
        
        config.ODEStepSize = 20;
        config.ODEType = 2;
        
        MSum = 0;
        numCells = 0;
        
        MSampled = 0;
        
        YMat = [];
        YMat_valid = [];
        
        validIdx = [];
        for i=1:length(prom.cellIdx)
            cellIdx = i;
            
            if (isempty(prom.TrReconstruction{i}) || prom.Valid(i)==0)
                continue;
            end
            
            
            numCells = numCells + 1;
            
            validIdx(numCells) = prom.cellIdx(i);
            %YMat(numCells, :) = prom.YFP{i}.Measurement;
        end
        
        
        cellIdx = 1:length(YFP.cells);
        
        for i=1:length(prom.cellIdx)
           YMat(i, :) = YFP.cells{prom.cellIdx(i)}.Measurement;%YFP.cells{i}.Measurement;
        end
        
        for i=1:length(validIdx)
           YMat_valid(i, :) = YFP.cells{validIdx(i)}.Measurement;
        end
        
        size(YMat_valid, 1) / size(YMat, 1)
        
%         figure(k);
%         subplot(rows, cols, u);
%         
%         plot(YMat', 'r'); hold on;
%         plot(YMat_valid', 'k');
%         
        
        MSum = MSum / numCells;
        
        muProteinData = mean(YMat);
        varProteinData = var(YMat);
        stdProteinData = sqrt(varProteinData);
        
        muProteinData_valid = mean(YMat_valid);
        varProteinData_valid = var(YMat_valid);
        stdProteinData_valid = sqrt(varProteinData_valid);
        
        figure(k);
        subplot(rows, cols, u);
        
        q = [];
        q(1, :) = muProteinData_valid - stdFactor*stdProteinData_valid;
        q(2, :) = 2*stdFactor*stdProteinData_valid;
        
        p = area(time/60, q'); hold on;
        set(p(1), 'FaceColor', 'none');
        set(p(1), 'LineStyle', 'none');
        set(p(2), 'FaceColor', 'g');
        set(p(2), 'LineStyle', 'none');
        set(p(2), 'FaceAlpha', 0.2);
        plot(time/60, muProteinData_valid, '-'); hold on;
        
        
        q = [];
        q(1, :) = muProteinData - stdFactor*stdProteinData;
        q(2, :) = 2*stdFactor*stdProteinData;
        
        p = area(time/60, q'); hold on;
        set(p(1), 'FaceColor', 'none');
        set(p(1), 'LineStyle', 'none');
        set(p(2), 'FaceColor', 'r');
        set(p(2), 'LineStyle', 'none');
        set(p(2), 'FaceAlpha', 0.2);
        plot(time/60, muProteinData, '-'); hold on;
        
        title(conditions{u}.Name, 'Interpreter', 'none');
        drawnow;
        
    end
    
end

