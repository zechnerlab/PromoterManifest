clear;
close all;

simulate=0;
addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');


if (simulate == 1)
    concVec = [100, 275, 690, 3000];
    conditions = GenerateConditions();
    %conditions = GetConditions(conditions, [50], concVec);
    conditions = conditions(20);
    promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
    %promoters = {'TKL2'};
    %promoters={'pSIP18_mut6'};
    
    stdFactor = 2.58;
    
    for i=1:length(conditions)
        for k=1:length(promoters)
            
            promName = promoters{k};
            conditionName = conditions{i}.Name;
            totName = sprintf('%s_%s', promName, conditionName);
            strs = split(conditionName, '_');
            conditionModelName = ['DM_50min_' strs{end}];
            conditionModelName = conditionName;
            %conditionModelName = [num2str)conditions{i}.Concentration 'nM_DM'];
            
            d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
            YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);
            res = load(['results/ModelInference_' promName '_' conditionModelName '.mat']);
            res.modelOpt.odeInfos = load('odeInfosFullMoments.mat');
            
            %subplot(1,length(testConditions), i);
            
            %cellIdx = setdiff(1:length(YFP.cells), res.cellIdx);
            %cellIdx = res.cellIdx;
            
            %cellIdx = 1:length(YFP.cells);
            
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
            
            q = [];
            stdMean = sqrt(uncertainties(1, :));
            meanMean = moments(1, :);
            meanVar = moments(2, :);
            stdVar = sqrt(uncertainties(2, :));
            
            time = linspace(0, max(YFP.cells{1}.MeasurementTime), 300);
            time = YFP.cells{1}.MeasurementTime;
            [Stats] = RunModel(res.modelOpt, time, d.TFInputParams);
            
            numMoments = length(res.modelOpt.odeInfos.infos.MomentSystem{1}.dM);
            numStates = length(res.modelOpt.Z);
            
            baseIdx = (0:numStates-1)*numMoments + numStates;
            
            mRNAMean = sum(Stats(baseIdx+1, :));
            proteinMean = sum(Stats(baseIdx+2, :));
            proteinVar = sum(Stats(baseIdx+7, :)) - proteinMean.^2;
            mRNAVar = sum(Stats(baseIdx+4, :)) - mRNAMean.^2;
            
            varYP = proteinVar;
            meanYP = proteinMean;
            
            %% plotting
            figure(i);
            subplot(2, length(promoters), k);
            q(1, :) = meanMean - stdFactor*stdMean;
            q(2, :) = 2*stdFactor*stdMean;
            p = area(YFP.cells{1}.MeasurementTime/60, q'); hold on;
            set(p(1), 'FaceColor', 'none');
            set(p(1), 'LineStyle', 'none');
            set(p(2), 'FaceColor', 'g');
            set(p(2), 'LineStyle', 'none');
            set(p(2), 'FaceAlpha', 0.2);
            plot(YFP.cells{1}.MeasurementTime/60, meanMean, '-');
            title(promoters{k});
            xlim([0, 155]);
            plot(time/60, meanYP, '-');
            
            
            subplot(2, length(promoters), length(promoters) + k);
            q(1, :) = meanVar - stdFactor*stdVar;
            q(2, :) = 2*stdFactor*stdVar;
            p = area(YFP.cells{1}.MeasurementTime/60, q'); hold on;
            set(p(1), 'FaceColor', 'none');
            set(p(1), 'LineStyle', 'none');
            set(p(2), 'FaceColor', 'g');
            set(p(2), 'LineStyle', 'none');
            set(p(2), 'FaceAlpha', 0.2);
            plot(YFP.cells{1}.MeasurementTime/60, meanVar, '-');
            title(promoters{k});
            xlim([0, 155]);
            plot(time/60, varYP, '-');
            
            
            MeanDistanceRel(i, k) = mean(abs(meanMean - proteinMean)./stdMean);
            VarDistanceRel(i, k) = mean(abs(meanVar - proteinVar)./stdVar);
            
            MeanDistanceAbs(i, k) = mean(abs(meanMean - proteinMean)./meanMean);
            VarDistancAbs(i, k) = mean(abs(meanVar - proteinVar)./meanVar);
           
            Protein(i, k) = max(meanMean);
            
        end
        
    end
    

    save results/FittingAccuracy.mat;
    
else
    
    load results/FittingAccuracy.mat
end


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




