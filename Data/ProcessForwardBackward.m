clear;
close all;

warning off;

addpath('../Common');
addpath('../../Common/Statistics/');
addpath('../../Common/StochChemKin/');
addpath('../../Common/Models/');
addpath('../../Common/ODEs/dopri');
addpath('../../Common/ODEs/');
addpath('../../Common/');
addpath('../../Common/SpecialFunctions/');
addpath('../../Data');
addpath('../../Data/MSN2');
addpath('../');


odeInfos = load('odeInfos.mat');

measurementNoiseResults = load('MeasurementNoiseHXK1_DM_50min_3uM.mat');

conditionsAll = GenerateConditions();
%conditions = conditions(21:26);


promoters = {'DCS2'};
%promoters = {'TKL2'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [690];
for u=1:length(cVec)
    conc = cVec(u);
  
  
    conditions = GetConditions(conditionsAll, [50], conc);
    modelCondition = GetConditions(conditionsAll, [50], conc);
    
    for i=1:length(conditions)
        for k=1:length(promoters)
            
            
            promName = promoters{k};
            conditionName = conditions{i}.Name;
            totName = sprintf('%s_%s', promName, conditionName);
            
            d = load(['../../Data/MSN2/' promName '/' totName '_MSN2.mat']);
            YFP = load(['../../Data/MSN2/' promName '/' totName '_YFP.mat']);
            results{k, i} = load(['ModelInference_' promName '_' modelCondition{1}.Name '.mat']);
            %modelCondition{1}.Name
            
            %cellIdx = randperm(length(YFP.cells));
            %YFP.cells = YFP.cells(cellIdx(1:min(length(cellIdx), 3000)));
            
            if (strcmp(modelCondition{1}.Name, conditionName))
                cellIdx = setdiff(1:length(YFP.cells), results{k, i}.cellIdx);
            else
                cellIdx = 1:length(YFP.cells);
            end
            Data.Conditions{1}.cells = YFP.cells(cellIdx);
            
            for j=1:length(Data.Conditions{1}.cells)
                Data.Conditions{1}.cells{j}.InputParams = d.TFInputParams;
                figure(99);
                plot(Data.Conditions{1}.cells{j}.MeasurementTime/60, Data.Conditions{1}.cells{j}.Measurement, 'o-'); hold on;
            end
            drawnow;
            hold off;
            
            
            model = results{k, i}.modelOpt;
            model.MeasurementSigma = measurementNoiseResults.modelOpt.MeasurementSigma;
            model.odeInfos = load('odeInfos.mat');
            model.odeInfosBackward = load('odeInfosBackward.mat');
            
            optParams = GetModelParameters(model, results{k, i}.config);
            
            [L, postStats, postTimes] = FilterForwardBackward(optParams, model, results{k, i}.config, Data);
            %model.c(3)
            
            pTime = postTimes{1};
            TMax = pTime(end);
            
            sGrid = linspace(0, TMax, 300);
            
            maxR2Vec = [];
            maxT2Vec = [];
            maxR1Vec = [];
            maxT1Vec = [];
            TVec = [];
            maxTrVec = [];
            trRateVec = [];
            trMat = [];
            P1Mat = [];
            P2Mat = [];
            
            Valid = repmat(true, 1, length(Data.Conditions{1}.cells));
            
            
            
            
            for l=1:length(Data.Conditions{1}.cells)
                
                postStat = postStats{l};
                postTime = postTimes{l};
                
                postStat = SampleCTMPPathGrid_mex(real(postStat), postTime, sGrid);
                postTime = sGrid;
                
                %             if (sum(sum(imag(postStat)))>0)
                %                Valid(l)= false;
                %                continue;
                %             end
                %
                
                if sum(sum(~isreal(postStat)))>0
                    Valid(l) = false;
                    continue;
                end
                
                
                
                trRate = model.Z*postStat(1:3, :);
                
                maskFunction = and(postTime>0*60,postTime<(conditions{i}.PulseParameters(2)*60 + 15*60));
                [maxTr, maxIdx] = max(trRate.*maskFunction);
                maxT = postTime(maxIdx);
                
                TVec(l) = maxT;
                maxTrVec(l) = maxTr;
                
                trMat(l, :) = trRate;%.*maskFunction;
                
                P1Mat(l, :) = postStat(2, :);
                P2Mat(l, :) = postStat(3, :);
                
                
            end
            
            colVec = 'rgbykmc';
            
            figure(1);
            plot(postTime/60, mean(trMat(Valid, :)), 'Color', colVec(1+mod(k, length(colVec)))); hold on;
            
            
            figure(2);
            p = plot(mean(TVec(Valid))/60, mean(maxTrVec(Valid)), 'o', 'Color', colVec(1+mod(k, length(colVec))));
            set(p, 'MarkerSize', 5+i*2);
            hold on; drawnow;
            
            
            Promoters{k}.MaxTrVec{i} = maxTrVec;
            Promoters{k}.TrMat{i} = trMat;
            Promoters{k}.P1Mat{i} = P1Mat;
            Promoters{k}.P2Mat{i} = P2Mat;
            Promoters{k}.MaxTrTimeVec{i} = TVec;
            Promoters{k}.MeanMaxTr(i) = mean(maxTrVec);
            Promoters{k}.MaskFunction = maskFunction;
            Promoters{k}.Model{i} = model;
            %Promoters{k}.
            Promoters{k}.Valid{i} = Valid;
            Promoters{k}.YFP{i} = Data;
            Promoters{k}.cellIdx = cellIdx;
            
            
            fprintf('Processed condition %s for promoter %s.\n', conditionName, promName);
            if (i==1)
                legendStr{end+1} = promName;
            end
        end
        
    end
    legend(legendStr);
    
    save(['Processed_' num2str(conditions{1}.Concentration) 'nM_DM.mat']);
    
end
