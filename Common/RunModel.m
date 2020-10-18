function [Stats] = RunModel(options, time, inputParams)



% M0 = options.M0;
% 
% 
% ODEOpts.RelTol = 0.001;
% ODEOpts.AbsTol = 0.001;
% 
% 
% options.Update = 0;
% options.InputParams = inputParams;
% 
% funHandleDop = @(t, y) options.odeInfos.odeHandle(t, y, options);
% 
% [~, Stats] = ode45(funHandleDop, time, M0, ODEOpts);
% 
% 
% Stats = Stats';


c(1) = options.Q(2, 1);
c(2) = options.Q(1, 2);
c(3) = options.Q(3, 2);
c(4) = options.Q(2, 3);
c(5:6) = options.Z(2:end);
c(7:9) = options.c(2:end);
modelIn.c = c;
modelIn.Update = 0;
modelIn.InputParams = inputParams;


M0 = options.odeInfos.infos.DefaultInitialConditions;

X0 = [1, 0, 0, options.InitialMeans, options.muMeas0, options.H(1)];

for k=1:length(M0)
    
    specIdx = options.odeInfos.infos.MomentSystem{1}.SpeciesIdx{k};
    M0(k) = prod(X0(specIdx));
    
end

%infos.DefaultInitialConditions = M0;
%infos.MomentSystem = moments;


%X0 = zeros(55, 1);
M0(12) = options.H(1);
M0(25) = options.muMeas0^2 + options.varMeas0;
M0(27) = options.H(1)^2 + options.H(2);
M0(22) = options.InitialMeans^2 + options.InitialVars;
M0(46) = options.H(1);
M0(48) = options.H(2) + options.H(1)^2;
M0(55) = options.H(2) + options.H(1)^2;
M0(47) = options.muMeas0*options.H(1);
M0(26) = options.muMeas0*options.H(1);
M0(11) = options.muMeas0;

opts = {};%{'RelTol', 0.01};
opts = {'RelTol', 0.001};
[~, x] = ode15s(@ODEFullMoments, time, M0, opts, modelIn);
Stats = x';
end