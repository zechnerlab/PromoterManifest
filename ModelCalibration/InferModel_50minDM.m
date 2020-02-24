clear all;
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
addpath('../Data');
addpath('../Data/MSN2');

%parpool(40);

odeInfos = load('odeInfosFullMoments.mat');

conditions = GenerateConditions();
promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3', 'DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
conditions = GetConditions(conditions, [50], [100, 275, 690, 3000]);
%conditions = fliplr(conditions);


AllResults = {};%cell(length(promoters), length(testConditions));

firstModel = cell(length(promoters), 1);


config = struct;
config.UnknownQIdx{1} = [1, 2];
config.UnknownQIdx{2} = [2, 1];
config.UnknownQIdx{3} = [2, 3];
config.UnknownQIdx{4} = [3, 2];
config.UnknownZIdx = [2,3];

config.UnknownCIdx = [2, 3];%[2, 3, 4];
config.CPriors.A = [20, 20];%[8, 8, 8];%][3, 5, 3];
config.CPriors.B = config.CPriors.A./[0.0013; 0.05];%config.CPriors.A./[0.0013; 0.25; 1.6667e-05];
config.UnknownHIdx = [2];
config.UnknownMIdx = [];

config.InferHyper = 0;
config.QPrior = 'Gamma';
% 

            
Ic = repmat(length(conditions):-1:1, length(promoters), 1);
Ip = repmat((1:length(promoters))', 1, length(conditions));
IdxVec = [Ip(:), Ic(:)];
conditionName = conditions{IdxVec(1, 2)}.Name;

Results = cell(size(IdxVec, 1), 1);


%promName = promoters{IdxVec(1, 1)};
%Results{1} = inferModel(promName, conditionName, config, 50000);
%firstModel{1} = Results{1}.modelOpt;   
%vfModel = firstModel{1};

%config.UnknownMIdx = [1];

parfor k=1:length(promoters)
    promName = promoters{IdxVec(k, 1)};
    Results{k} = inferModel(promName, conditionName, config, 20000);
    firstModel{k} = Results{k}.modelOpt;   
end


config.UnknownCIdx = [];
config.CPriors.A = [];
config.CPriors.B = [];%config.CPriors.A./[0.0013; 0.25; 1.6667e-05];
config.UnknownHIdx = [];
config.UnknownMIdx = [];

parfor u=length(promoters)+1:size(IdxVec, 1)
    
    promName = promoters{IdxVec(u, 1)};
    conditionName = conditions{IdxVec(u, 2)}.Name;
    
    Results{u} = inferModel(promName, conditionName, config, 20000, firstModel{IdxVec(u, 1)});
    
    fprintf('Finished run %d (%d)!\n', u, size(IdxVec, 1));
end
      

for u=1:length(Results)
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
    conditionName = conditions{IdxVec(u, 2)}.Name;
    save(['results/ModelInference_' promName '_' conditionName '.mat'], 'model', 'modelOpt', 'chain', 'targetData', 'cellIdx', 'config', 'startParams', 'K', 'validIdx', 'promName', 'conditionName');
end

SaveInferenceStats;



