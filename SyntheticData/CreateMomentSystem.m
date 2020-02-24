clear all;
close all;

addpath('../Common');
addpath('../Common/ODEs/dopri/');
warning off;


[Pre, Post, c, X0] = CreateTwoStageExpressionSystemZ();
X0(1) = 0;
X0(2) = 0;
%X0(3) = 0;


modelIn.Pre = Pre;
modelIn.Post = Post;
modelIn.c = c;
modelIn.X0 = X0;
%modelIn.InputParams.inputTransitionIdx = [1, 2];
modelIn.ObservedSpeciesIdx = 2;

%[odeHandle, infos] = GenerateMomentsInput(modelIn);
%save odeInfosLN.mat odeHandle infos;

load('odeInfos.mat');

promName = 'SIP18';
resultsInference = load(['results/ModelInference_' promName '_DM_50min_275nM.mat']);

M = 100;
resultsInference.modelOpt.InputODEInfos.odeHandle = odeHandle;
resultsInference.modelOpt.InputODEInfos.infos = infos;
resultsInference.modelOpt.MeasurementSigma = 800;

d = load(['../Data/MSN2/' promName '/' promName '_DM_50min_275nM_MSN2.mat']);
YFP = load(['../Data/MSN2/' promName '/' promName '_DM_50min_275nM_YFP.mat']);

resultsInference.modelOpt.InputParams = d.TFInputParams;
model = resultsInference.modelOpt;
model.NumBins = 8;

filterLength = 59;

tGrid = linspace(0, max(YFP.cells{1}.MeasurementTime), 100);

numCells = length(YFP.cells);

meanRna = zeros(numCells, length(tGrid));
meanProtein = zeros(numCells, length(tGrid));
meanZ = zeros(numCells, length(tGrid));
varRna = zeros(numCells, length(tGrid));
varProtein = zeros(numCells, length(tGrid));
varZ = zeros(numCells, length(tGrid));

MeanTauS = zeros(numCells, 1);
MeanTauA = zeros(numCells, 1);
PActivated = zeros(numCells, 1);

tic

parfor j=1:3
    
    
    cellStruct = YFP.cells{j};
    
    particles = InitializeParticleDistribution(model, M);
    particles = ReconstructCell(model, particles, cellStruct, filterLength);
    
    [Mu, Sigma, pActivated, P2, P3, meanTauS, meanTauA] = ProcessCellReconstruction(model, particles, tGrid);
    
    meanRna(j, :) = Mu(1, :);
    meanProtein(j, :) = Mu(2, :);
    meanZ(j, :) = Mu(3, :);
    
    varRna(j, :) = Sigma(1, :);
    varProtein(j, :) = Sigma(2, :);
    varZ(j, :) = Sigma(3, :);
    
    MeanTauS(j) = meanTauS;
    MeanTauA(j) = meanTauA;
    PActivated(j) = pActivated;
    
    fprintf('Finished cell %d (tauS=%f, tauA=%f, pA=%f)\n', j, meanTauS/60, meanTauA/60, pActivated);
    
    
    plot(tGrid, meanProtein(j, :), tGrid, meanProtein(j, :) - sqrt(varProtein(j, :)), tGrid, meanProtein(j, :) + sqrt(varProtein(j, :)), cellStruct.MeasurementTime, cellStruct.Measurement, 'o'); drawnow;
end

toc


