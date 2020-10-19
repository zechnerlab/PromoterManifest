clear;
close all;

addpath('../Common');
addpath('../StateReconstruction/');

generate = 1;

[Pre, Post, c, X0] = CreateTwoStageExpressionSystemZ();

modelIn.Pre = Pre;
modelIn.Post = Post;
modelIn.c = c;
modelIn.X0 = X0;
modelIn.ObservedSpeciesIdx = 2;

if (generate == 1)
    [odeHandle, infos] = GenerateMomentsInputLN(modelIn);
    save odeInfosLN.mat odeHandle infos;
else
    load odeInfosLN.mat;
end


