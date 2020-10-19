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


concVec = [100, 275, 690, 3000];
conditions = GenerateConditions();
conditions = GetConditions(conditions, [50], concVec);

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

for k=1:length(promoters)
    for i=1:length(conditions)
        
        
        promName = promoters{k};
        conditionName = conditions{i}.Name;
        totName = sprintf('%s_%s', promName, conditionName);
        strs = split(conditionName, '_');
        conditionModelName = ['DM_50min_' strs{end}];
        
        res = load(['results/ModelInference_' promName '_' conditionModelName '.mat']);
        
        Z(i, :) = res.modelOpt.Z;
        
    end
    
    maxZ = max(max(Z));
    
    promStats{k}.maxZ = maxZ;
end

save 'results/promStats.mat';

