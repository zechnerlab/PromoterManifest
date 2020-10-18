function [particle] = SimulatePromoter(model, particlePrior, T)


M = length(particlePrior);

startTime = particlePrior{1}.t(end);
particle = particlePrior;

if (nargin<3)
    T = max(model.InputParams.inputTimes);
end


syncTimes = linspace(startTime, T, 100);

for k=1:M
    
    
    
    
    % initialize model parameters for particular particle
    [Pre, Post, c, X0] = CreatePromoterModel(particle{k}.Model.Q);
    
    inputParams = model.InputParams;
    inputParams.inputLevels = inputParams.inputLevels*c(1);
    inputParams.inputTimes = [0, inputParams.inputTimes];

    X0(:) = 0;
    X0(particlePrior{k}.stateIdx(end)) = 1;
    
    [X, tTmp, r] = SimulateSSA_TV(X0, c, Pre, Post, 60000, T, 1, inputParams, [inputParams.inputTimes, T, T], startTime);
    
    t = tTmp(2:end);
    jIdx = find(r~=0);
    pT = t(jIdx);
    pR = r(jIdx);
    stateIdx = sum(repmat((1:length(X0))', 1, size(jIdx, 2)).*X(:, 1+jIdx));
    
    particle{k}.t = [particle{k}.t, pT];
    %stateIdx2 = sum(repmat((1:length(X0))', 1, size(X, 2)-1).*X(:, 2:end));
    particle{k}.stateIdx = [particle{k}.stateIdx, stateIdx];
    particle{k}.Z = [particle{k}.Z, particle{k}.Model.Z(stateIdx)];
    particle{k}.r = [particle{k}.r, pR];
    
    particle{k}.Z = [particle{k}.Z, particle{k}.Z(end)];
    particle{k}.t = [particle{k}.t, tTmp(end)];
    particle{k}.r = [particle{k}.r, 0];
    particle{k}.stateIdx = [particle{k}.stateIdx, particle{k}.stateIdx(end)];
end
%stairs(t, stateIdx); hold on;
%plot(inputParams.inputTimes, inputParams.inputLevels(1:end-1)/max(inputParams.inputLevels));
end