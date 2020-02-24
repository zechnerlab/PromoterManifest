function particle = InitializeParticleDistribution(model, config, M, posteriorSamples)

if (nargin<4)
   posteriorSamples = [];
   priorSampling = 0;
else
   priorSampling = 1; 
end

M0 = model.InputODEInfos.infos.DefaultInitialConditions;

model.InitialMeans = 0;
model.InitialVars = 0;
means = [model.InitialMeans; model.muMeas0; model.H(1)];
vars = [model.InitialVars; model.varMeas0; model.H(2)];

M0 = SetInitialInputMoments(model.InputODEInfos.infos, M0, means, vars);

parIdx = randperm(size(posteriorSamples, 2));

particle = cell(M, 1);

for i=1:M
    particle{i}.t = 0;
    particle{i}.stateIdx = (rand>model.P0(1)) + 1;
    particle{i}.Z = model.Z(particle{i}.stateIdx);
    particle{i}.r = [];
    particle{i}.PosteriorStats = M0;
    particle{i}.PosteriorTime = 0;
    particle{i}.L = 0;
    particle{i}.Idx = i;
    if (priorSampling == 1)
        particle{i}.Model = SetModelKineticParamsOnly(model, config, posteriorSamples(:, parIdx(i)));
    else
        particle{i}.Model = model;
    end
    
end

end