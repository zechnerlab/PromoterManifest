%% PlotTemporalReconstructions.m 
%Example script to demonstrate how the inferred transcription and promoter 
%state dynamics can be plotted over time. 

clear;
close all;

addpath('../StateReconstruction');
addpath('../StateReconstruction/results');
addpath('../Common/');
addpath('../Data/Msn2');

%Set the vector of all promoter names
promoters = {'DCS2', 'HXK1', 'SIP18', 'ALD3','DDR2', 'RTN2', 'TKL2', 'pSIP18_mut6', 'pSIP18_mut21'};

%Specify which promoter you want to analyze and set corresponding index.
selPromoter = 'SIP18';
k=find(strcmp(promoters, selPromoter));

%set this one to 1 if you want to plot responding cells only, or zero
%otherwise.
respondersOnly = 0;

%Generate vector of all 30 conditions. 
conditions = GenerateConditions();

%Set pulse length (10, 20, ..., 50) and inducer concentration (100, 275, 690,
%3000) and retrieve the corresponding condition from the 'conditions'
%structure. Note that only the single pulse experiments can be obtained using the function
%'GetConditions'. Repeated pulse experiments have indices 21-30, those have to be selected
%manually.
pulseLength = 50;
concentration = 690;
condition = GetConditions(conditions, pulseLength, concentration);

%specify which run to show.
resultsName = 'results_1';
results = load(['../StateReconstruction/' resultsName '/StateReconstruction_' condition{1}.Name '.mat'], 'Promoters');


if (respondersOnly==1)
    selectedCells = and(results.Promoters{k}.PActivated>0.99, results.Promoters{k}.Valid);
else
    selectedCells = logical(results.Promoters{k}.Valid);
end

meanTr = mean(results.Promoters{k}.MeanTr(selectedCells, :), 1);
trRec = results.Promoters{k}.TrReconstruction(selectedCells);

% calculate average promoter state probabilities
P2 = mean(results.Promoters{k}.P2(selectedCells, :));
P3 = mean(results.Promoters{k}.P3(selectedCells, :));
P1 = 1 - P2 - P3;

semFactor = 2.58; %specify how many std-devs the error bars should correspond to


q(1, :) = meanTr - semFactor*std(results.Promoters{k}.MeanTr(selectedCells, :), 1)/sqrt(sum(selectedCells));
q(2, :) = 2*semFactor*std(results.Promoters{k}.MeanTr(selectedCells, :), 1)/sqrt(sum(selectedCells));

t = results.Promoters{k}.tGrid; %Set time vector

subplot(1,2,1);
plot(t/60, P1, t/60, P2, t/60, P3);
xlabel('Time in min');
ylabel('Promoter state probability');
legend('State 0', 'State 1', 'State 2');

subplot(1,2,2);
p = area(t/60, q'); hold on;
set(p(1), 'FaceColor', 'none');
set(p(1), 'EdgeColor', 'none');
set(p(2), 'FaceAlpha', 0.1, 'FaceColor', 'r', 'EdgeColor', 'r', 'EdgeAlpha', 0.2);
plot(t/60, meanTr, 'Color',  'r');
xlabel('Time in min');
ylabel('Transcription rate in 1/s');
xlim([0, 100]);
p = title([promoters{k} ' (' condition{1}.Name ')']);
set(p, 'Interpreter', 'none');


