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


concVec = [100, 275, 690, 3000];
conditions = GenerateConditions();

promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

plotCondition = conditions{30};


promName = 'SIP18';
conditionName = plotCondition.Name;
totName = sprintf('%s_%s', promName, conditionName);

d = load(['../Data/MSN2/' promName '/' totName '_MSN2.mat']);
YFP = load(['../Data/MSN2/' promName '/' totName '_YFP.mat']);

M = 200;
cellIdx = randperm(length(YFP.cells));

for l=1:M
    
    YMat(l, :) = YFP.cells{cellIdx(l)}.Measurement;
    
end


p = plot(YFP.cells{1}.MeasurementTime/60, YMat); hold on;
set(p, 'Color', [0.0, 0.2, 0.4, 0.4]);
box off;





