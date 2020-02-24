clear;
close all;


addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');


concVec = [100, 275, 690, 3000];
conditions = GenerateConditions();
%conditions = GetConditions(conditions, [50], concVec);
conditions = conditions(20);


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'TKL2'};
%promoters={'pSIP18_mut6'};

validation = 1;

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
        
        if (validation == 1) %show cells *not* used for fitting
            cellIdx = setdiff(1:length(YFP.cells), res.cellIdx);
        else %show cells used for fitting
            cellIdx = res.cellIdx;
        end
        
        %cellIdx = 1:length(YFP.cells);
        
        YMat = [];
        for l=1:length(cellIdx)
            
            YMat(l, :) = YFP.cells{cellIdx(l)}.Measurement;

        end
        
        for l=1:size(YMat, 2)
           data{l} = YMat(:, l); 
        end
        
        %[moments, uncertainties] = BootstrapMomentsStd(data, 1000);
        
        q = [];
        stdY = std(YMat);
        meanY = mean(YMat);
        stdFactor = 1;
        q(1, :) = meanY - stdFactor*stdY;
        q(2, :) = 2*stdFactor*stdY;
        figure(i);
        subplot(1,length(promoters), k);
        p = area(YFP.cells{1}.MeasurementTime/60, q'); hold on;
        set(p(1), 'FaceColor', 'none');
        set(p(1), 'LineStyle', 'none');
        set(p(2), 'FaceColor', 'g');
        set(p(2), 'LineStyle', 'none');
        set(p(2), 'FaceAlpha', 0.2);
        plot(YFP.cells{1}.MeasurementTime/60, meanY, '-');
        title(promoters{k});
        xlim([0, 155]);
        %p = plot(YFP.cells{1}.MeasurementTime/60, YMat(randi(length(cellIdx), 20, 1), :), '-g'); hold on;
        %set(p, 'Color', [0, 0, 0, 0.2]);
        
        time = linspace(0, max(YFP.cells{1}.MeasurementTime), 300);
        %res.modelOpt.H(2) = 0.01;
        [Stats] = RunModel(res.modelOpt, time, d.TFInputParams);
        

        if (i==1)
           res.modelOpt.c;
        end
        
        res.modelOpt.Z%model.Z
        
        numMoments = length(res.modelOpt.odeInfos.infos.MomentSystem{1}.dM);
        numStates = length(res.modelOpt.Z);
        
        baseIdx = (0:numStates-1)*numMoments + numStates;
        
        mRNAMean = sum(Stats(baseIdx+1, :));
        proteinMean = sum(Stats(baseIdx+2, :));
        proteinVar = sum(Stats(baseIdx+7, :)) - proteinMean.^2;
        mRNAVar = sum(Stats(baseIdx+4, :)) - mRNAMean.^2;
        
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
        
%         [P, R] = RunModelSSA(res.modelOpt, time, d.TFInputParams, 100);
%         plot(time/60, mean(P));
        
        condNames{i} = conditionName;
    end
    
end

legend(condNames);




