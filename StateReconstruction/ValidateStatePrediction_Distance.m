clear;
close all;

simulate = 1;
addpath('../Common/');
addpath('../Data/MSN2/');

if simulate == 1
    
    promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
    selProm = promoters;%{'SIP18'};
    %selProm = {'SIP18'};
    
    legendStr = {};
    cVec = [100, 275, 690, 3000];
    plVec = [10, 20, 30, 40, 50];
    
    conditions = GenerateConditions();
    conditions = conditions(20);
    
    
    rows = 5;
    cols = 9;
    
    stdFactor = 1;
    
    for u=1:length(conditions)
        
        res = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
        
        for j=1:length(selProm)
            
            promName = selProm{j};
            k = find(strcmp(promoters, promName));
            
            totName = sprintf('%s_%s', promName, conditions{u}.Name);
            YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
            
            
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
            
            %M0 = SetInitialInputMoments(model.InputODEInfos.infos, M0, means, vars);
            
            config.ODEStepSize = 10;
            config.ODEType = 2;
            
            MSum = 0;
            numCells = 0;
            
            MSampled = 0;
            
            YMat = [];
            PMat = [];
            
            for i=1:length(prom.cellIdx)
                cellIdx = i;
                
                if (isempty(prom.TrReconstruction{i}) || prom.Valid(i)==0)
                    continue;
                end
                
                if (sum(isnan(prom.meanProtein(i, :)))>0)
                    continue;
                end
                
                
                if (sum(abs(imag(prom.meanProtein(i, :)))>0)>0)
                    continue;
                end
                
                means = [model.InitialMeans; model.muMeas0; prom.meanZ(cellIdx, end)];
                vars = [model.InitialVars; model.varMeas0; 0];
                
                M0 = model.InputODEInfos.infos.DefaultInitialConditions;
                M0 = SetInitialInputMoments(model.InputODEInfos.infos, M0, means, vars);
                
                trProfile.t = prom.TrReconstruction{cellIdx}.t;
                trProfile.Value = prom.TrReconstruction{cellIdx}.Z;
                
                
                output = SolveGeneModel(prom.tGrid, model, config, trProfile, M0);
                %plot(output.t, output.M(2, :), ...
                %    prom.tGrid, prom.meanProtein(cellIdx, :), '-b'); hold on; drawnow;
                
                MSampled = SampleCTMPPathGrid_mex(output.M, output.t, time);
                MSum = MSum + MSampled;
                
                
                
                numCells = numCells + 1;
                %YMat(numCells, :) = prom.YFP{i}.Measurement;
                PMat(numCells, :) = MSampled(2, :);
            end
            
            
            MSum = MSum / numCells;
            
            muProtein = MSum(2, :);
            varProtein = MSum(7, :) - muProtein.^2;
            stdProtein = sqrt(varProtein);
            
            
            cellIdx = 1:length(YFP.cells);
            
            YMat = [];
            for l=1:length(cellIdx)
                
                YMat(l, :) = YFP.cells{cellIdx(l)}.Measurement;
                
            end
            
            data = cell(0, 0);
            
            for l=1:size(YMat, 2)
                data{l} = YMat(:, l);
            end
            
            [moments, uncertainties, q5, q95] = BootstrapMomentsStd(data, 1000);
            
            muProteinData = moments(1, :);
            varProteinData = moments(2, :);
            stdProteinData = sqrt(varProteinData);
            muProteinDataStd = sqrt(uncertainties(1, :));
            varProteinDataStd = sqrt(uncertainties(2, :));
            
            MeanDistanceRel(u, j) = mean(abs(muProtein - muProteinData)./muProteinDataStd);
            VarDistanceRel(u, j) = mean(abs(varProtein - varProteinData)./varProteinDataStd);
            
            MeanDistanceAbs(u, j) = mean(abs(muProtein - muProteinData)./muProteinData);
            VarDistanceAbs(u, j) = mean(abs(varProtein - varProteinData)./varProteinData);
            
            Protein(u, j) = max(muProteinData);
            
            figure(u);
            subplot(rows, cols, k);
            
            
            q = [];
            q(1, :) = muProteinData - stdFactor*stdProteinData;
            q(2, :) = 2*stdFactor*stdProteinData;
            
            p = area(time/60, q'); hold on;
            set(p(1), 'FaceColor', 'none');
            set(p(1), 'LineStyle', 'none');
            set(p(2), 'FaceColor', 'g');
            set(p(2), 'LineStyle', 'none');
            set(p(2), 'FaceAlpha', 0.2);
            plot(time/60, muProteinData, '-g'); hold on;
            
            q = [];
            q(1, :) = muProtein - stdFactor*stdProtein;
            q(2, :) = 2*stdFactor*stdProtein;
            
            p = area(time/60, q'); hold on;
            set(p(1), 'FaceColor', 'none');
            set(p(1), 'LineStyle', 'none');
            set(p(2), 'FaceColor', 'r');
            set(p(2), 'LineStyle', 'none');
            set(p(2), 'FaceAlpha', 0.2);
            plot(time/60, muProtein, '-r'); hold on;
            
            
            
            
            title(conditions{u}.Name, 'Interpreter', 'none');
            drawnow;
            
        end
        
    end
    
    %save results/FittingAccuracy_SMC.mat;
else
    load results/FittingAccuracy_SMC.mat;
end

figure;
numBins = 30;
subplot(2,3,1);
[h, b] = hist(log10(MeanDistanceRel(:)), numBins);
h = h/sum((b(2) - b(1))*h);
area(b, h, 'FaceAlpha', 0.3, 'FaceColor', 'b'); hold on;

[h, b] = hist(log10(VarDistanceRel(:)), numBins);
h = h/sum((b(2) - b(1))*h);
area(b, h, 'FaceAlpha', 0.3, 'FaceColor', 'r');

subplot(2,3,2);
plot(log10(Protein), log10(MeanDistanceRel), '.b');

subplot(2,3,3);
plot(log10(Protein), log10(VarDistanceRel), '.b');

subplot(2,3,4);
bar(median(MeanDistanceRel));

subplot(2,3,5);
bar(median(VarDistanceRel));
