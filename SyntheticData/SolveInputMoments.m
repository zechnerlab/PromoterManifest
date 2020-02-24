function [Mu, Sigma, tGrid] = SolveInputMoments(model, t, Z)

    model.TranscriptionRateSequence.inputTimes = t(2:end);
    model.TranscriptionRateSequence.inputLevels = Z(1:end);
    
    M0 = model.InputODEInfos.infos.DefaultInitialConditions;
    means = [model.InitialMeans; model.muMeas0; model.H(1)];
    vars = [model.InitialVars; model.varMeas0; model.H(2)];
    
    M0 = SetInitialInputMoments(model.InputODEInfos.infos, M0, means, vars);
    
    tGrid = linspace(0, max(t), 200);
    model.Update = 0;
    [~, M] = ode45(model.InputODEInfos.odeHandle, tGrid, M0, {'RelTol', 1e-2, 'AbsTol', 1e-4}, model);
    M = M';
    
    Mu = M(1:3, :);
    Sigma = M(4:9, :);
    

    
    %subplot(1,3,3);
    %plot([0, model.TranscriptionRateSequence.inputTimes], model.TranscriptionRateSequence.inputLevels);
end