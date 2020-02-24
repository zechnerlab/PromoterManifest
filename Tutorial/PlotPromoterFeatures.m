%% PlotPromoterFeatures.m 
%Example script to demonstrate how arbitrary promoter features can be
%plotted against each other for all single pulse conditions. 

clear;
close all;

addpath('../StateReconstruction');
addpath('../StateReconstruction/results');
addpath('../Common/');
addpath('../Data/Msn2');


%load promoter data obtained from single-cell trajectory inference.
data = load('../StateReconstruction/results/PromoterFeatures.mat');
dataBase = data.dataB; %dataBase is the core structure that stores all promoter features.

%extract promoter names, concentrations, number of pulses, 
[promoters, concentrations, numPulses, durations, intervals] = GetDataParameters(dataBase);

%The 'concentrations' vector contains the concentrations of the inducer
%molecule. 'Msn2level' contains the corresponding Msn2
%induction levels that are reached upon treatment with the respective
%inducer concentration.
Msn2level = [0.25, 0.5, 75, 1];

%Generate different colors used to indicate different Msn2 induction
%levels.
colVec = GetDefaultColors();

%Contains the names of all transcriptional features that were calculated
%for each of the promoters.
features = dataBase.featureNames;
numFeatures = length(features);

%String used for constructing the plot legend.
legendStr = {};


%Sets the name of the promoter that you want to analyze. Possible promoters
%are 'DCS2', 'HXK1', 'SIP18', 'ALD3', 'DDR2', 'TKL2', 'RTN2',
%'pSIP18_mut6', 'pSIP18_mut21'. Note that 'pSIP18_mut6' corresponds to
%'pSIP18 A4' and 'pSIP18_mut21' corresponds 'pSIP18 D6' in the main text.
promName = 'DDR2';
l = find(strcmp(promoters, promName));

%Select the two promoter features that you want to plot against each other.
%Possible features are 
%    'pActivate' (percentage of responders)
%    'meanTauS' (mean time to activate)
%    'meanTauA' (mean of time active)
%    'numResponders' (number of responding cells)
%    'numCells' (number of valid cells that were processed)
%    'timeToMaxTr' (time until the average transcription rate reaches a
%    maximum)
%    'maxTr' (maximum of the average transcription rate)
%    'meanTrOutput' (average transcriptional output)
%    'meanTrOutputR' (average transcriptional output considering only
%    responding cells)
feature1 = 'meanTrOutputR';
feature2 = 'meanTauS';


%specify how many standard errors the error bars should indicate.
semFactor = 2.58;

p = zeros(length(concentrations), 1);

for k=1:length(concentrations)
    
    %read out the values of the two features for a given concentration and
    %all pulse lengths
    [data1] = GetDataMatrix(dataBase, promoters{l}, concentrations(k), [1], [], [-1], feature1);
    [data2] = GetDataMatrix(dataBase, promoters{l}, concentrations(k), [1], [], [-1], feature2);
    
    
    %read out standard error of each feature. Standard errors are available
    %for all features except the percentage of responders.
    if (~strcmp(feature1, 'pActivate'))
        sem1 = [feature1 'Std'];
        [data1sem] = GetDataMatrix(dataBase, promoters{l}, concentrations(k), [1], [], [-1], sem1);
    else
       data1sem = zeros(size(data1)); 
    end
    if (~strcmp(feature2, 'pActivate'))
        sem2 = [feature2 'Std'];
        [data2sem] = GetDataMatrix(dataBase, promoters{l}, concentrations(k), [1], [], [-1], sem2);  
    else
       data2sem = zeros(size(data2)); 
    end

    %plot the feature1 against feature2
    p(k) = plot(data1, data2, '-'); hold on;
    set(p(k), 'Color', colVec(k, :));
    Plot2DErrorbar(data1, data2, data1sem, data2sem, [colVec(k, :), 0.8], semFactor);
    
    %mark the individual points as colored discs of varying size. Larger
    %pulse-lengths have larger disc size.
    for j=1:length(data1)
        p2 = plot(data1(j), data2(j), 'o');
        set(p2, 'MarkerFaceColor', colVec(k, :), 'MarkerEdgeColor', colVec(k, :)); 
        set(p2, 'MarkerSize', 2 + j*1);
    end
    
    %Set labels according to feature1 and feature2. The title of the plot
    %is just the promoter name
    title(promoters{l});
    xlabel(feature1);
    ylabel(feature2);
 
    legendStr{k} = num2str(concentrations(k));
end
legend(p, legendStr);

