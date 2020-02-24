function model = CreateDefaultModel(numStates, config, TFInputParams, muMeas0, varMeas0)


model.Q = 0.01*ones(numStates, numStates);
model.Q = model.Q - diag(sum(model.Q, 1));
model.Z = 0.01*rand(numStates, 1);
%k1 = 0.002;
%k2 = 0.008;

k1 = 0.001;
k2 = 0.001;

model.Q = [0,  k2 0.001
            k1, 0 0.001
            0.001   0.001  0];
model.Q = model.Q - diag(sum(model.Q, 1));
%         
% model.Q = [-0.02, 0.01 0.01
%             0.01, -0.02 0.01
%             0.01   0.01  -0.02];


%model.Z(1) = 0.01;% = 0.1*zeros(numStates, 1);
model.Z = [0.0001, 0.01, 0.00001];
model.c = [0; 0.0013; 0.05; 1.6667e-05];

muExt = 0.4;
sigmaExt_squared = 0.01;
model.H(1) = muExt;
model.H(2) = sigmaExt_squared;
model.MeasurementSigma = 60;
val = 10e-3;
if (strcmp(config.QPrior, 'LaplaceAll'))
    model.alpha = val*rand(length(config.UnknownQIdx), 1);
    model.beta = val*rand(length(config.UnknownZIdx), 1);
elseif (strcmp(config.QPrior, 'LaplaceGroup'))
    model.alpha = val*rand(numStates, 1);
    model.beta = val*rand(length(config.UnknownZIdx), 1);
elseif (strcmp(config.QPrior, 'GammaAll'))
    model.alpha = 1;
    model.beta = 10;
else
    model.alpha = val;
    model.beta  = val;
end

%model.alpha = 10e-3*ones(length(config.UnknownQIdx), 1);
%model.beta = 10e-3*ones(length(config.UnknownZIdx), 1);
model.InputParams = TFInputParams;
model.InputParams.inputLevels = TFInputParams.inputLevels + eps;

odeInfos = load('odeInfos.mat');
X0 = odeInfos.infos.DefaultInitialConditions;
X0(1:numStates) = 1/numStates;
X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], 3, muExt*1/numStates); %set E[Z]
X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [3 3], (muExt^2+sigmaExt_squared)*1/numStates); %set E[ZZ]
X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [2], muMeas0*1/numStates);
X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [2 2], (muMeas0^2 + varMeas0)*1/numStates);
X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [2 3], muMeas0*muExt*1/numStates);

model.muMeas0 = muMeas0;
model.varMeas0 = varMeas0;

model.M0 = X0;
model.P0 = X0(1:numStates);
model.NumBins = 35;
model.odeInfos = odeInfos;
model.StatusOutput = 0;
end