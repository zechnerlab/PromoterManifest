function model = CreateDefaultModel(numStates, config, TFInputParams, muMeas0, varMeas0)


model.Q = 0.01*ones(numStates, numStates);
model.Q = model.Q - diag(sum(model.Q, 1));
model.Z = 0.01*rand(numStates, 1);


k = 0.01;

        
model.Q = [ 0         k       0
            k         0       k
            0         k      0];
model.Q = model.Q;

model.Q = model.Q - diag(sum(model.Q, 1));

model.Z = [0, 0.1, 0.2]; 
model.c = [0; 0.0013; 0.025; 1.6667e-05;];


muExt = 1;
sigmaExt_squared = 0.2;
model.H(1) = muExt;
model.H(2) = sigmaExt_squared;
model.MeasurementSigma = 0;
% for ln model
%model.MeasurementSigma = 0.1;


model.InputParams = TFInputParams;
model.InputParams.inputLevels = TFInputParams.inputLevels + eps;

odeInfos = load('odeInfosFullMoments.mat');
X0 = odeInfos.infos.DefaultInitialConditions;

model.NumSpecies = sum(odeInfos.infos.MomentSystem{1}.Order==1);

X0(1:numStates) = [1, 0, 0];

model.InitialMeans = [0;];
model.InitialVars = [0;];


model.muMeas0 = muMeas0;
model.varMeas0 = varMeas0;

%model.M0 = X0;
model.P0 = X0(1:numStates);
model.NumBins = 25;
model.odeInfos = odeInfos;
model.StatusOutput = 0;
end