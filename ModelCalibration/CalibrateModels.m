%This file calibrates models for each promoter and conditions separately.
%Only the promoter paramters (switching rates, transcription rates) are
%allwed to vary between conditions. The remaining parameters are estimated
%from the 50min-pulse 100% Msn2 condition. The results are stored in the
%folder 'results'. Once the calibration is finished, the function
%'SaaveInferenceStats.m' is executed, which calculates the maximal transcription
%rate for from the 50min pulse conditions for each promoter respectively.
%This is used later during state reconstruction to distinguish responding from 
%non-responding cells. Note that the script takes a few hours to complete on 
%a 40-core shared memory system.

clear;
close all;

addpath('../Common');
addpath('../Data');
addpath('../Data/MSN2');


odeInfos = load('odeInfosFullMoments.mat');

conditions = GenerateConditions();
promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3', 'DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

concVec = [100, 275, 690, 3000];

conditions50min3uM = GetConditions(conditions, [50], [3000]);
conditions50min100_690nM = GetConditions(conditions, [50], [100, 275, 690]);

conditionsRest = GetOtherConditions(conditions, [conditions50min3uM, conditions50min100_690nM]);


Results = cell(length(conditions)*length(promoters), 1);

%% Perform calibration of complete model using the 50min 3uM pulse condition
%Both the promoter parameters (theta) and the parameters associated with
%mRNA degradation, translation and protein degr
firstModel = cell(length(promoters), length(concVec));

%Number of samples used for MCMC sampling
numSamples = 20000;
 
config = struct;
config.UnknownQIdx{1} = [1, 2];
config.UnknownQIdx{2} = [2, 1];
config.UnknownQIdx{3} = [2, 3];
config.UnknownQIdx{4} = [3, 2];
config.UnknownZIdx = [2,3];

config.UnknownCIdx = [2, 3];
config.CPriors.A = [20, 20];
config.CPriors.B = config.CPriors.A./[0.0013, 0.05];
config.UnknownHIdx = [2];
config.UnknownMIdx = [];

config.QPrior = 'Gamma'; 


usedConditions = conditions50min3uM;
Ic = repmat(1:1:length(usedConditions), length(promoters), 1);
Ip = repmat((1:length(promoters))', 1, length(usedConditions));
IdxVec = [Ip(:), Ic(:)];
conditionName = usedConditions{IdxVec(1, 2)}.Name;

concIdx = find(concVec == 3000); %for 3uM

parfor k=1:size(IdxVec, 1)
    promName = promoters{IdxVec(k, 1)};
    Results{k} = inferModel(promName, conditionName, config, numSamples);
    Results{k}.condition = usedConditions{IdxVec(k, 2)};
    
    firstModel{k, concIdx} = Results{k}.modelOpt;   %store fitted model for the subsequent inference runs
    
    fprintf('Finished run %d (%d)!\n', k, size(IdxVec, 1));
end


%% Perform calibration for all 50min conditions with 25%, 50% and 75% Msn2.
%In this case not all parameters are re-estimated, but only the promoter
%switching rates + transcription rates (mRNA / protein decays + translation
%are used from the inital run).

config.UnknownCIdx = [];
config.CPriors.A = [];
config.CPriors.B = [];
config.UnknownHIdx = [];
config.UnknownMIdx = [];

idxOffset = size(IdxVec, 1);

usedConditions = conditions50min100_690nM;
Ic = repmat(1:1:length(usedConditions), length(promoters), 1);
Ip = repmat((1:length(promoters))', 1, length(usedConditions));
IdxVec = [Ip(:), Ic(:)];
conditionName = usedConditions{IdxVec(1, 2)}.Name;

firstModelIdx = concIdx; %use parameter estimates from initial fit (50min, 100%Msn2)
numSamples = 20000;

parfor u=1:size(IdxVec, 1)
    
    promName = promoters{IdxVec(u, 1)};
    conditionName = usedConditions{IdxVec(u, 2)}.Name;
    
    Results{idxOffset + u} = inferModel(promName, conditionName, config, numSamples, firstModel{IdxVec(u, 1), firstModelIdx});
    Results{idxOffset + u}.condition = usedConditions{IdxVec(u, 2)};
     
    fprintf('Finished run %d (%d)!\n', u, size(IdxVec, 1));
end

%store into firstModelStructure for subsequent runs
for u=1:size(IdxVec, 1)
   concIdx = find(usedConditions{IdxVec(u, 2)}.Concentration == concVec);
   firstModel{IdxVec(u, 1), concIdx} = Results{idxOffset + u}.modelOpt;
end


%% Perform calibration for all remaining conditions 
%As before, only the switching and transcription rates are reestimated.

config.UnknownCIdx = [];
config.CPriors.A = [];
config.CPriors.B = [];
config.UnknownHIdx = [];
config.UnknownMIdx = [];

idxOffset = idxOffset + size(IdxVec, 1);

usedConditions = conditionsRest;
Ic = repmat(1:1:length(usedConditions), length(promoters), 1);
Ip = repmat((1:length(promoters))', 1, length(usedConditions));
IdxVec = [Ip(:), Ic(:)];
conditionName = usedConditions{IdxVec(1, 2)}.Name;

parfor u=1:size(IdxVec, 1)
    
    promName = promoters{IdxVec(u, 1)};
    conditionName = usedConditions{IdxVec(u, 2)}.Name;
    
    concIdx = find(usedConditions{IdxVec(u, 2)}.Concentration == concVec);
    
    Results{idxOffset + u} = inferModel(promName, conditionName, config, numSamples, firstModel{IdxVec(u, 1), concIdx});
    Results{idxOffset + u}.condition = usedConditions{IdxVec(u, 2)};
    
    fprintf('Finished run %d (%d)!\n', u, size(IdxVec, 1));
end
      
      

for u=1:length(Results)
    if (~isempty(Results{u}))
        model = Results{u}.model;
        modelOpt = Results{u}.modelOpt;
        chain = Results{u}.chain;
        targetData = Results{u}.targetData;
        cellIdx = Results{u}.cellIdx;
        config = Results{u}.config;
        K = Results{u}.K;
        validIdx = Results{u}.validIdx;
        promName = Results{u}.promName;
        startParams = Results{u}.startParams;
        conditionName = Results{u}.condition.Name;
        save(['results/ModelInference_' promName '_' conditionName '.mat'], 'model', 'modelOpt', 'chain', 'targetData', 'cellIdx', 'config', 'startParams', 'K', 'validIdx', 'promName', 'conditionName');
    end
end

SaveInferenceStats;


