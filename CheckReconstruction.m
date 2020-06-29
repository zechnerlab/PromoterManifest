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


odeInfos = load('odeInfos.mat');

%measurementNoiseResults = load('MeasurementNoiseHXK1_DM_50min_3uM.mat');

conditionsAll = GenerateConditions();
%conditions = conditions(21:26);


promoters = {'HXK1', 'DCS2', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};
%promoters = {'DCS2'};

counter = 0;
maxCond = 8;

legendStr = {};
cVec = [100, 275, 690, 3000];

conditions = conditionsAll;

colVec = 'rgbkc';

condIdx = 15;

u = condIdx;

if (conditions{u}.PulseParameters(1) == 1) %% single pulse
    conc = conditions{u}.Concentration;
    modelCondition =  GetConditions(conditionsAll, [50], conc);
else %% repeated pulses
    modelCondition =  GetConditions(conditionsAll, [50], 690);
end

results = load(['results/StateReconstruction_' conditions{u}.Name '.mat'], 'Promoters');

cIdx = find(cVec == conc);
plIdx = round(conditions{u}.PulseParameters(2)/10);

pIdx = 3;

k = pIdx;


numC = 50;

tGrid = linspace(0, max(results.Promoters{k}.YFP{1}.MeasurementTime), size(results.Promoters{k}.meanRna, 2));
responders = results.Promoters{k}.PActivated>0.5';
%valid = and(sum(imag(results.Promoters{k}.meanRna)>0, 2)==0, sum(real(results.Promoters{k}.meanRna<0)==0, 2));
valid = results.Promoters{k}.Valid;

for l=1:numC*numC
    

% 
% 
% figure(1);
% subplot(numC, numC, l);
% plot(tGrid, real(results.Promoters{k}.meanRna(l, :)'));
% 
% figure(2);
% subplot(numC, numC, l);
% p = plot(tGrid, real(results.Promoters{k}.MeanTr(l, :)));
% set(p, 'LineWidth', 2);
% ylim([0, max(results.Promoters{k}.Model.Z)]);
% 
% 
% figure(3);
% subplot(numC, numC, l);
p = plot(results.Promoters{k}.YFP{l}.MeasurementTime, results.Promoters{k}.YFP{l}.Measurement); hold on;

if (valid(l)==1)
    if (responders(l)==1)
        set(p, 'Color', 'r');
    else
        set(p, 'Color', 'k');
    end
else
    set(p, 'Color', 'g');
end



%drawnow;


end

