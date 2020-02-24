clear all;
close all;

addpath('../../../Common/CurveFitting');
addpath('../../Common/');

promoters = EnumeratePromoters('.');

for l=1:length(promoters)
    promoterName = promoters{l};
    
    addpath(promoterName);
    
    
    conditions = GenerateConditions();
    
    meanYFP = [];
    varYFP = [];
    
    for k=1:length(conditions)
        
        expName = [promoterName '_' conditions{k}.Name];
        
        data = load([expName '.mat']);
        
        
        startIdx = 1;
        Time = data.time(startIdx:end);
        YFP = data.YFP(:, startIdx:end);
        CFP = data.CFP(:, startIdx:end);
        MSN2 = data.MSN2_RFP(:, startIdx:end);
        
        T = max(Time);
        numCells = size(YFP, 1);
        
        [gridTF, gridT] = MSN2_trace_generator(conditions{k}.Concentration,...
            [0, T], conditions{k}.PulseParameters);
        
        maxTF = max(gridTF);
        TFInputParams.inputLevels = gridTF/1000;% / maxTF;
        TFInputParams.inputTimes = gridT(2:end)*60;
        TFInputParams.inputRateIndex = 1;
        TFInputParams.Func = @(t)MSN2_trace_generator_oneVal(conditions{k}.Concentration, [t/60], conditions{k}.PulseParameters)/maxTF;
        
        %     subplot(1,2,1);
        %     plot(TFInputParams.inputTimes / 60, gridTF(1:end-1), 'r');hold on;
        %     plot(Time, mean(MSN2), 'b');
        %
        
        
        
        maturationIdx = 1;
        Time = Time(maturationIdx:end);
        Time = Time - min(Time);
        meanYFP = mean(YFP(:, 1:maturationIdx), 2);
        YFP = YFP(:, maturationIdx:end);
        %YFP = YFP - repmat(meanYFP, 1, size(YFP, 2));
        
        %validIdx = sum(YFP==0, 2)==0;
        %YFP = YFP(validIdx, :);
        numCells = size(YFP, 1);
        
        offset = 1;
        cells = {};
        currIdx = 1;
        for i=1:numCells
            Measurement = YFP(i, offset+1:end);
            if (sum(Measurement)>50)
                cells{currIdx}.Measurement = YFP(i, offset+1:end);
                cells{currIdx}.MeasurementTime = Time(offset+1:end)*60;
                cells{currIdx}.Simulated = 0;
                currIdx = currIdx + 1;
            end
            
%             if (i<20 && k==18)
%                 %subplot(1,2,2)
%                 p = plot(cells{i}.MeasurementTime, cells{i}.Measurement, 'o--k'); hold on;
%                 set(p, 'MarkerFaceColor', 'r');
%                 set(p, 'MarkerEdgeColor', 'k');
%             end

        end
        
        Promoter{l}.Conditions{k}.Cells = cells;
        Promoter{l}.Conditions{k}.name = conditions{k}.Name;
        
        save([promoterName '/' promoterName '_' conditions{k}.Name '_YFP.mat'], 'cells');
        save([promoterName '/' promoterName '_' conditions{k}.Name '_MSN2.mat'], 'TFInputParams');
        
    end
    
    fprintf('Processed promoter %s...\n', promoterName); 
    
end