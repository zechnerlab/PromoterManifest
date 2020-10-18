clear;
close all;

simulate = 0;
validation = 1; %set to one if comparison should be shown against data that has *NOT* been used for fitting

addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');


resultNames = { 'results_1'
                'results_2'
                'results_3'
                'results_4'
                'results_5'};

stdFactor = 1;

if (simulate == 1)
    concVec = [100, 275, 690, 3000];
    conditions = GenerateConditions();
    %conditions = GetConditions(conditions, [50], concVec);
    
    plotCondition = conditions(20); %plot 50min 100% Msn2 condition as an example
    promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

    
    for j=1:length(resultNames)
        resultName = resultNames{j};
        
        for i=1:length(conditions)
            for k=1:length(promoters)
                
                promName = promoters{k};
                conditionName = conditions{i}.Name;
                totName = sprintf('%s_%s', promName, conditionName);
                strs = split(conditionName, '_');
                
                conditionModelName = conditionName;
                
                d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
                YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
                res = load([resultName '/ModelInference_' promName '_' conditionModelName '.mat']);
                res.modelOpt.odeInfos = load('odeInfosFullMoments.mat');
                
                
                if (validation == 1) %show cells *not* used for fitting
                    cellIdx = setdiff(1:length(YFP.cells), res.cellIdx);
                else %show cells used for fitting
                    cellIdx = res.cellIdx;
                end
                
                YMat = [];
                for l=1:length(cellIdx)
                    
                    YMat(l, :) = YFP.cells{cellIdx(l)}.Measurement;
                    
                end
                
                YMatAll = [];
                
                for l=1:length(YFP.cells)
                    
                    YMatAll(l, :) = YFP.cells{l}.Measurement;
                    
                end
                
                data = cell(0, 0);
                
                for l=1:size(YMat, 2)
                    data{l} = YMat(:, l);
                end
                
                [moments, uncertainties, q5, q95] = BootstrapMomentsStd(data, 1000);
                
                q = [];
                stdMean = sqrt(uncertainties(1, :));
                meanMean = moments(1, :);
                meanVar = moments(2, :);
                stdVar = sqrt(uncertainties(2, :));
                
                time = YFP.cells{1}.MeasurementTime;
                [Stats] = RunModel(res.modelOpt, time, d.TFInputParams);
                
                numMoments = length(res.modelOpt.odeInfos.infos.MomentSystem{1}.dM);
                numStates = length(res.modelOpt.Z);
                
                %Read out the corresponding moments (only protein needed but
                %RNA is computed as well);
                mRNAMean = Stats(4, :);
                proteinMean = Stats(5, :);
                proteinVar = Stats(25, :) - proteinMean.^2;
                mRNAVar = Stats(22, :) - mRNAMean.^2;
                
                varYP = proteinVar;
                meanYP = proteinMean;
                
                %% plotting
                if (strcmp(conditionName, plotCondition{1}.Name))
                    %figure(j);
                    subplot(length(resultNames), length(promoters), (j-1)*length(promoters) + k);
                    q = [];
                    stdY = std(YMat);
                    meanY = mean(YMat);
                    
                    q(1, :) = meanY - stdFactor*stdY;
                    q(2, :) = 2*stdFactor*stdY;

                    p = area(YFP.cells{1}.MeasurementTime/60, q'); hold on;
                    set(p(1), 'FaceColor', 'none');
                    set(p(1), 'LineStyle', 'none');
                    set(p(2), 'FaceColor', 'g');
                    set(p(2), 'LineStyle', 'none');
                    set(p(2), 'FaceAlpha', 0.2);
                    plot(YFP.cells{1}.MeasurementTime/60, meanY, '-');
                    title(promoters{k});
                    xlim([0, 155]);
                    
                    
                    q = [];
                    stdYP = sqrt(proteinVar);
                    meanYP = proteinMean;
                    q(1, :) = meanYP - stdFactor*stdYP;
                    q(2, :) = 2*stdFactor*stdYP;
                    
                    p = area(time/60, q'); hold on;
                    set(p(1), 'FaceColor', 'none');
                    set(p(1), 'LineStyle', 'none');
                    set(p(2), 'FaceColor', 'r');
                    set(p(2), 'LineStyle', 'none');
                    set(p(2), 'FaceAlpha', 0.2);
                    plot(time/60, meanYP, '-');
                    
                    
                    drawnow;
                end
                
                %Calculate relative distance between model and data
                %moments.
                MeanDistanceRel(i, k, j) = mean(abs(meanMean - proteinMean)./stdMean);
                VarDistanceRel(i, k, j) = mean(abs(meanVar - proteinVar)./stdVar);
                
                %store maximum of the average protein copy number calculated
                %over all cells. This is used to test whether there is a
                %correlation between the model error and the protein abundance.
                Protein(i, k, j) = max(mean(YMatAll));
                
                fprintf('Finished condition %d, promoter %d, run %d\n', i, k, j);
            end
            
        end
    end
    
    
    save resultsAvrg/FittingAccuracy.mat;
    
else
    
    load resultsAvrg/FittingAccuracy.mat
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




