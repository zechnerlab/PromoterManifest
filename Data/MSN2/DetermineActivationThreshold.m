clear all;
close all;

addpath('../../../Common/CurveFitting');
addpath('../../Common/');

promoters = EnumeratePromoters('.');

%for l=1:length(promoters)
    promoterName = 'HXK1';
    
    addpath(promoterName);
    
    
    conditions = GenerateConditions();
    conditions = conditions(20);
    
    
    meanYFP = [];
    varYFP = [];
    
    for k=1:length(conditions)
        
        expName = [promoterName '_' conditions{k}.Name '_size'];
        
        data = load([expName '.mat']);
        
        
        startIdx = 1;
        Time = data.time(startIdx:end);
        YFP = data.YFP_molecules(:, startIdx:end);
        CFP = data.CFP_molecules(:, startIdx:end);
        MSN2 = data.MSN2_RFP(:, startIdx:end);
        
        subplot(1,2,1);
        plot(Time, mean(YFP), 'o-');
      
        subplot(1,2,2);
        plot(Time, mean(MSN2));
        
    end
    
    fprintf('Processed promoter %s...\n', promoterName); 
    
%end