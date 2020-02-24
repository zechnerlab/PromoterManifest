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
%conditions = conditions(21:26);


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'DCS2'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;

colVec = 'rgbkcmy';

dataB = struct;

numCond = 20;
for u=1:20%numCond%length(conditions)
    
    results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    
    if (conditions{u}.PulseParameters(1) == 1)
        maxTimeWindow = conditions{u}.PulseParameters(2) + 10;
    else
        maxTimeWindow = conditions{u}.PulseParameters(1)*(conditions{u}.PulseParameters(2) + conditions{u}.PulseParameters(3))+10;
    end
    

    for k=1:length(promoters)
        [features] = ExtractPromoterFeatures(results.Promoters{k}, 1, maxTimeWindow*60);
        dataB = StoreFeaturesToDatabase(dataB, features, promoters{k}, conditions{u});
        
        fprintf('Stored promoter %s, condition %s\n', promoters{k}, conditions{u}.Name);
    end

end

save results/PromoterFeatures.mat dataB;



