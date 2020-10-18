clear all;
close all;

addpath('../Common');
addpath('../Common/ODEs/dopri/');
addpath('../StateReconstruction/');
warning off;

generate = 1;

[Pre, Post, c, X0] = CreateSystemForFullMoments();

%X0(1) = 0;
%X0(2) = 0;
%X0(3) = 0;


modelIn.Pre = Pre;
modelIn.Post = Post;
modelIn.c = c;
modelIn.X0 = X0;
%modelIn.InputParams.inputTransitionIdx = [1, 2];
modelIn.ObservedSpeciesIdx = 4;

if (generate == 1)
    [odeHandle, infos] = GenerateFullMoments(modelIn);
    save odeInfosFullMoments.mat odeHandle infos;
else
    load odeInfosFullMoments.mat;
end

modelIn.Update = 0;
modelIn.InputParams.inputLevels = [1, 1];
modelIn.InputParams.inputTimes = [0, 200];
[t, x] = ode15s(odeHandle, [0, 150*60], infos.DefaultInitialConditions, {}, modelIn);
x = x';
subplot(1,3,1);
plot(t, x(1, :), t, x(2, :), t, x(3, :));

subplot(1,3,2);
plot(t, x(4, :), t, x(5, :));


