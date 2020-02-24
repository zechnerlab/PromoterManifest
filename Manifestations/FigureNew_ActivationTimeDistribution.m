clear all;
close all;

addpath('results/');


promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();

colVec = 'rgbkcmy';


k=find(strcmp(promoters, 'SIP18'));

for u=1:20
    
    results = load(['../StateReconstruction/results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');
    t = results.Promoters{k}.tGrid;
    
    concIdx = find(conditions{u}.Concentration == cVec);
    plIdx = floor(conditions{u}.PulseParameters(2)/10);
    
    responders = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
    
    tauS = results.Promoters{k}.MeanTauS(responders);
    tauA = results.Promoters{k}.MeanTauA(responders);
    
    [h, b] = hist(tauA, 10);
    db = b(2) - b(1);
    h = h/(db*sum(h));
    
    figure(99);
    subplot(length(cVec), 1+length(plVec), (concIdx-1)*(1+length(plVec)) + 1);
    scatter(tauS/60, tauA/60, '.'); hold on;
    %stairs(b/60, h); hold on;
    title([promoters{k} '(' num2str(conditions{u}.Concentration) 'nM)']);
    
    
    
    drawnow;
    
    fprintf('Processed promoter %s, condition %s\n', promoters{k}, conditions{u}.Name);
    
end





