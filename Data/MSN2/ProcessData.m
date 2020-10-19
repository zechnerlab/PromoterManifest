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
        
        expName = [promoterName '_' conditions{k}.Name '_size'];
        
        data = load([expName '.mat']);
        
        startIdx = 1;
        Time = data.time(startIdx:end);
        
        %divide copy numbers by cell size (gives ~concentration) and then multiply by the average
        %cell size per time point. This should reduce errors due to segmentation.
        YFP = data.YFP_molecules(:, startIdx:end)./data.cell_size_pixels(:, startIdx:end) .* mean(data.cell_size_pixels(:, startIdx:end));
        
        CFP = data.CFP_molecules(:, startIdx:end);
        MSN2 = data.MSN2_RFP(:, startIdx:end);
        
        T = max(Time);
        numCells = size(YFP, 1);
        
        [gridTF, gridT] = MSN2_trace_generator(conditions{k}.Concentration,...
            [0, T], conditions{k}.PulseParameters);
        
        maxTF = max(gridTF);
        TFInputParams.inputLevels = gridTF/1000; %divide by some value to make closer to 1.
        %This rescaling is only for plotting purposes and is irrelevant for
        %the analysis since the input is multiplied by a free parameter
        %gamma, which is estimated from the data.
        TFInputParams.inputTimes = gridT(2:end)*60;
        TFInputParams.inputRateIndex = 1;
        
        %also function handle to the Msn2 trace generator is stored.
        %However this is obselete and not used anymore.
        TFInputParams.Func = @(t)MSN2_trace_generator_oneVal(conditions{k}.Concentration, [t/60], conditions{k}.PulseParameters)/maxTF;
        
        %shift trajectories by 12.5min, roughtly corresponding to
        %maturation delay.
        maturationIdx = 6;
        Time = Time(maturationIdx:end);
        Time = Time - min(Time);
        meanYFP = mean(YFP(:, 1:maturationIdx), 2);
        YFP = YFP(:, maturationIdx:end);
        CFP = CFP(:, maturationIdx:end);
        
        numCells = size(YFP, 1);
        
        offset = 0;
        cells = {};
        currIdx = 1;
        for i=1:numCells
            Measurement = YFP(i, offset+1:end);
            TF = MSN2(i, offset+1:end);
            
            if (sum(Measurement)>=50) %remove cells which had essentially no signal.
                cells{currIdx}.Measurement = YFP(i, offset+1:end);
                cells{currIdx}.MeasurementTime = Time(offset+1:end)*60;
                cells{currIdx}.Simulated = 0;
                currIdx = currIdx + 1;
            end
        end
        
        Promoter{l}.Conditions{k}.Cells = cells;
        Promoter{l}.Conditions{k}.name = conditions{k}.Name;
        
        save([promoterName '/' promoterName '_' conditions{k}.Name '_YFP.mat'], 'cells');
        save([promoterName '/' promoterName '_' conditions{k}.Name '_MSN2.mat'], 'TFInputParams');
        
    end
    
    fprintf('Processed promoter %s...\n', promoterName);
    
end