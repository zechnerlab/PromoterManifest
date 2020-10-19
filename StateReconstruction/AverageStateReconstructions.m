clear;
close all;

addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');


conditionsAll = GenerateConditions();

resultsName{1} = 'results_1';
resultsName{2} = 'results_2';
resultsName{3} = 'results_3';
resultsName{4} = 'results_4';
resultsName{5} = 'results_5';


for i=1:length(resultsName)
    
    load([resultsName{i} '/PromoterFeatures.mat']);
    
    [promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataB);
    
    dataAll(:, i) = dataB.Data(:);
    DataCube(:, :, i) = dataB.Data;
end

dataMean = nanmean(dataAll, 2);
dataSEM = nanstd(dataAll, [], 2)./sqrt(sum(~isnan(dataAll), 2)); 

dataMean = reshape(dataMean, size(dataB.Data));
dataSEM = reshape(dataSEM, size(dataB.Data));
dataB.Data = dataMean;
dataB.DataSEM = dataSEM;
dataB.DataCube = DataCube;
mkdir('resultsAvrg');
save resultsAvrg/PromoterFeatures.mat dataB;