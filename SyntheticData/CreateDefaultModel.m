function model = CreateDefaultModel(numStates, TFInputParams, muMeas0, varMeas0)


model.Q = 0.01*ones(numStates, numStates);
model.Q = model.Q - diag(sum(model.Q, 1));
model.Z = 0.01*rand(numStates, 1);
%k1 = 0.002;
%k2 = 0.008;

k1 = 0.01;
k2 = 0.0002;
k = 0.01;

model.Q = k1*ones(numStates, numStates);
model.Q = model.Q - diag(sum(model.Q, 1));


        
model.Q = [ 0         k       0
            k         0       k
            0         k      0];
model.Q = model.Q + 0*0.0001;

model.Q = model.Q - diag(sum(model.Q, 1));

%model.Z(1) = 0.01;% = 0.1*zeros(numStates, 1);
model.Z = [0, 0.1, 0.2]; %0.01*ones(1, numStates);%[0.01, 0.01, 0.01];
model.c = [0; 0.0013; 0.025; 1.6667e-05;];
%model.c = [0; 0.0013; 0.25; 0.0012; 1.6667e-05;];

muExt = 1;
sigmaExt_squared = 0.2;
model.H(1) = muExt;
model.H(2) = sigmaExt_squared;
model.MeasurementSigma = 0;
% for ln model
%model.MeasurementSigma = 0.1;

 
%model.alpha = 10e-3*ones(length(config.UnknownQIdx), 1);
%model.beta = 10e-3*ones(length(config.UnknownZIdx), 1);
model.InputParams = TFInputParams;
model.InputParams.inputLevels = TFInputParams.inputLevels + eps;

odeInfos = load('odeInfosFullMoments.mat');
X0 = odeInfos.infos.DefaultInitialConditions;

model.NumSpecies = sum(odeInfos.infos.MomentSystem{1}.Order==1);

X0(1:numStates) = [1, 0, 0];

model.InitialMeans = [0;];
model.InitialVars = [0;];

means = [model.InitialMeans; muMeas0; muExt];
vars = [model.InitialVars; varMeas0; sigmaExt_squared];


X0 = SetInitialMoments(odeInfos.infos, X0, numStates, means, vars);

model.muMeas0 = muMeas0;
model.varMeas0 = varMeas0;

model.M0 = X0;
model.P0 = X0(1:numStates);
model.NumBins = 25;
model.odeInfos = odeInfos;
model.StatusOutput = 0;
end