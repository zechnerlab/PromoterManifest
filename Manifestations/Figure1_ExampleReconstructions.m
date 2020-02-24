clear all;
close all;

addpath('results/');
addpath('../Data/Msn2');

load ../StateReconstruction/results/PromoterFeatures.mat;

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'DCS2'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];
plVec = [10, 20, 30, 40, 50];

conditions = GenerateConditions();

colVec = 'rgbkcmy';

k=find(strcmp(promoters, 'SIP18'));




plotCondition = conditions{19};
results = load(['../StateReconstruction/results/StateReconstruction_' plotCondition.Name '.mat'], 'Promoters');

%results.Promoters{k}.
%cellIdx = [6, 12]-1; %19 85
cellIdx = [2, 7];
numReconstructions = length(cellIdx);

area([0, results.Promoters{1}.YFP{1}.InputParams.inputTimes]/60, results.Promoters{1}.YFP{1}.InputParams.inputLevels);

figure;
for i=1:numReconstructions
    j = cellIdx(i);
    
    if (results.Promoters{k}.Valid(j)==1)
    subplot(numReconstructions, 3, (i-1)*3 + 1);
    
    
    
    plot(results.Promoters{k}.YFP{j}.MeasurementTime/60, results.Promoters{k}.YFP{j}.Measurement, '.-'); hold on;
    
    subplot(numReconstructions, 3, (i-1)*3 + 2);
    
    stairs(results.Promoters{k}.tGrid/60, results.Promoters{k}.MeanTr(j, :), '-'); hold on;
    %stairs(results.Promoters{k}.TrReconstruction{j}.t, results.Promoters{k}.TrReconstruction{j}.Z); 
    %stairs(results.Promoters{k}.tGrid/60, results.Promoters{k}.MeanTr(j, :) - sqrt(results.Promoters{k}.VarTr(j, :)) , '-');
    %stairs(results.Promoters{k}.tGrid/60, results.Promoters{k}.MeanTr(j, :) + sqrt(results.Promoters{k}.VarTr(j, :)) , '-');
    ylims = ylim;
    ylims(1) = 0;
    ylim(ylims);
    
    subplot(numReconstructions, 3, (i-1)*3 + 3);
    
    
    
    dT = (results.Promoters{k}.tGrid(2) - results.Promoters{k}.tGrid(1))/(2*60);
    
    t = [results.Promoters{k}.tGrid; results.Promoters{k}.tGrid];
    t = t(:);
    
    p2 = [results.Promoters{k}.P2(j, :); results.Promoters{k}.P2(j, :)];
    p2 = p2(:);
    
    p3 = [results.Promoters{k}.P3(j, :); results.Promoters{k}.P3(j, :)];
    p3 = p3(:);
    
    a = area(t(2:end)/60, p2(1:end-1)); hold on;
    set(a, 'FaceAlpha', 0.3, 'FaceColor', 'r');
    %set(a, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'BarWidth', 1);
    %stairs(results.Promoters{k}.tGrid/60, results.Promoters{k}.P2(j, :), 'Color', [0.5, 0.5, 0.5]);
    
    a = area(t(2:end)/60, p3(1:end-1)); hold on;
    set(a, 'FaceAlpha', 0.3, 'FaceColor', 'b');
    %set(a, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'BarWidth', 1);
    %stairs(results.Promoters{k}.tGrid/60, results.Promoters{k}.P3(j, :), 'Color', [0.5, 0.5, 0.5]);
    ylim([0, 1.5]);
    end
end






