clear;
close all;

addpath('../Common/');
addpath('../Data/MSN2/');

simulate = 0;
validation = 1;

resultsNames = { 'results_1'
    'results_2'
    'results_3'
    'results_4'
    'results_5'};

if simulate == 1
    
    promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
    selProm = promoters;
    
    legendStr = {};
    cVec = [100, 275, 690, 3000];
    plVec = [10, 20, 30, 40, 50];
    
    conditions = GenerateConditions();
    plotCondition = conditions(20); %plot 50min 100% Msn2 condition as an example

    stdFactor = 1;
    
    
    for o=1:length(resultsNames)
        resultsName = resultsNames{o};
        
        for u=1:length(conditions)
            
            res = load(['../StateReconstruction/' resultsName '/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
            
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
                
                %calculate average over the non-central moments
                MSum = MSum / numCells;
                
                %convert to central moments;
                muProtein = MSum(2, :);
                varProtein = MSum(7, :) - muProtein.^2;
                stdProtein = sqrt(varProtein);
                
                if (validation == 1) %show cells *not* used for fitting
                    cellIdx = setdiff(1:length(YFP.cells), prom.cellIdx);
                else %show cells used for fitting
                    cellIdx = prom.cellIdx;
                end
                
                YMat = [];
                for l=1:length(cellIdx)
                    YMat(l, :) = YFP.cells{cellIdx(l)}.Measurement;
                end
                
                %also collect all trajectories in a structure. This is needed
                %only to compute the total average concentration.
                YMatAll = [];
                for l=1:length(YFP.cells)
                    YMatAll(l, :) = YFP.cells{l}.Measurement;
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
                
                MeanDistanceRel(u, j, o) = mean(abs(muProtein - muProteinData)./muProteinDataStd);
                VarDistanceRel(u, j, o) = mean(abs(varProtein - varProteinData)./varProteinDataStd);
                
                Protein(u, j, o) = max(mean(YMatAll));
                
                %% plotting
                if (strcmp(conditions{u}.Name, plotCondition{1}.Name))
                    
                    
                    subplot(length(resultsNames), length(selProm), (o-1)*length(selProm) + j);
                    
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
                fprintf('Finished condition %d, promoter %d, run %d\n', u, j, o);
            end
            
        end
    end
    
    save resultsAvrg/FittingAccuracy_SMC.mat;
else
    load resultsAvrg/FittingAccuracy_SMC.mat;
end

MeanDistanceRel = mean(MeanDistanceRel, 3);
VarDistanceRel = mean(VarDistanceRel, 3);

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
plot(log10(Protein(:, :, 1)), log10(MeanDistanceRel), '.b');

subplot(2,3,3);
plot(log10(Protein(:, :, 1)), log10(VarDistanceRel), '.b');

subplot(2,3,4);
bar(median(MeanDistanceRel));

subplot(2,3,5);
bar(median(VarDistanceRel));

