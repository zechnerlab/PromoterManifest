clear;
close all;

warning off;

addpath('../Common');
addpath('../Common/Statistics/');
addpath('../Common/StochChemKin/');
addpath('../Common/Models/');
addpath('../Common/ODEs/dopri');
addpath('../Common/ODEs/');
addpath('../Common/');
addpath('../Common/SpecialFunctions/');
addpath('../Data');
addpath('../Data/MSN2');
addpath('../');


conditionsAll = GenerateConditions();

% resultsName{1} = 'results 08-Oct-2020 20:42:34';
% resultsName{2} = 'results 09-Oct-2020 06:46:30';
% resultsName{3} = 'results 09-Oct-2020 16:57:03';
% resultsName{4} = 'results 09-Oct-2020 22:25:49';
% resultsName{5} = 'results 10-Oct-2020 03:33:11';

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