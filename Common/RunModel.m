function [Stats] = RunModel(options, time, inputParams)



M0 = options.M0;


ODEOpts.RelTol = 0.001;
ODEOpts.AbsTol = 0.001;


options.Update = 0;
options.InputParams = inputParams;

funHandleDop = @(t, y) options.odeInfos.odeHandle(t, y, options);

[~, Stats] = ode45(funHandleDop, time, M0, ODEOpts);

%[Stats] = ode2(funHandleDop, time, M0);

Stats = Stats';
end