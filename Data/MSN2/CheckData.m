clear all;
close all;

promoters = EnumeratePromoters('.');
conditions = GenerateConditions;

promoterIdx = 4;
promoter = promoters{promoterIdx};

load([promoter '/' promoter '_' conditions{12}.Name '_YFP.mat']);
load([promoter '/' promoter '_' conditions{12}.Name '_MSN2.mat']);
for i=1:length(cells)
   
    subplot(1,2,1);
    plot(TFInputParams.inputTimes, TFInputParams.inputLevels(1:end-1), 'r-.');
    
    
    subplot(1,2,2);
    plot(cells{i}.MeasurementTime, cells{i}.Measurement, 'o--k'); hold on;
    
end