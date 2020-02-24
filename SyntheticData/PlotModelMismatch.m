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

sVec = [0, 0.1, 0.25];

perturbationType = 2;

for u=1:length(sVec)
    
    perturbationStrength = sVec(u);
    
    fileName = sprintf('FakeDataReconstruction_%f_%d.mat', perturbationStrength, perturbationType);
    res = load(fileName);
    
    for k=1:length(res.Promoters)
        
        prom = res.Promoters{k};
        modelTmp = CreateDefaultModel(3, prom.Model.InputParams, 0, 0);
        
        %plot true model moments
        prom.Model.M0 = modelTmp.M0;
        prom.Model.odeInfos = modelTmp.odeInfos;
        Stats = RunModel(prom.Model, prom.tGrid, prom.Model.InputParams);
        
        numMoments = length(prom.Model.odeInfos.infos.MomentSystem{1}.dM);
        numStates = length(prom.Model.Z);
        
        baseIdx = (0:numStates-1)*numMoments + numStates;
        
        mRNAMean = sum(Stats(baseIdx+1, :));
        proteinMean = sum(Stats(baseIdx+2, :));
        proteinVar = sum(Stats(baseIdx+7, :)) - proteinMean.^2;
        mRNAVar = sum(Stats(baseIdx+4, :)) - mRNAMean.^2;
        
        subplot(length(sVec),length(res.Promoters), (u-1)*length(res.Promoters) + k);
        plot(prom.tGrid, proteinMean, 'r', prom.tGrid, proteinMean - sqrt(proteinVar), 'r',...
            prom.tGrid, proteinMean + sqrt(proteinVar), 'r'); hold on;
        
        %plot perturbed model moments
        model = prom.Models{1};
        model.M0 = modelTmp.M0;
        model.odeInfos = modelTmp.odeInfos;
        Stats = RunModel(model, prom.tGrid, model.InputParams);
        
        mRNAMean = sum(Stats(baseIdx+1, :));
        proteinMean = sum(Stats(baseIdx+2, :));
        proteinVar = sum(Stats(baseIdx+7, :)) - proteinMean.^2;
        mRNAVar = sum(Stats(baseIdx+4, :)) - mRNAMean.^2;
        
        
        plot(prom.tGrid, proteinMean, 'b', prom.tGrid, proteinMean - sqrt(proteinVar), 'b',...
            prom.tGrid, proteinMean + sqrt(proteinVar), 'b'); hold on; 
        drawnow;
    end
end
