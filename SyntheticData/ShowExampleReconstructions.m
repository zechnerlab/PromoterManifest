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
sVec = 0;
numShownCells = 20;
cols = 6;
rows = ceil(numShownCells / cols);

for u=1:length(sVec)
    
    perturbationStrength = sVec(u);
    
    fileName = sprintf('FakeDataReconstruction_%f.mat', perturbationStrength);
    res = load(fileName);
    
    for k=1:length(res.Promoters)
        
        figure;
        prom = res.Promoters{k};
        
        numCells = size(prom.MeanTauA, 1);
        
        for l=1:numShownCells
            
            meanTr = prom.MeanTr(l, :);
            stdTr = sqrt(prom.VarTr(l, :));
            
            subplot(rows, cols, l);
            p = plot(prom.tGrid, meanTr, 'r'); hold on;
            plot(prom.tGrid, meanTr - stdTr, 'r');
            plot(prom.tGrid, meanTr + stdTr, 'r');
            
            stairs(prom.YFP{l}.t, prom.YFP{l}.Z, 'b');
        end
        
       
    end
end
