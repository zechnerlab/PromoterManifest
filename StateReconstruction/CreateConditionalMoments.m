clear all;
close all;

addpath('../Common');
addpath('../Common/ODEs/dopri/');
addpath('../StateReconstruction/');
warning off;

generate = 1;

[Pre, Post, c, X0] = CreateTwoStageExpressionSystemZ();



modelIn.Pre = Pre;
modelIn.Post = Post;
modelIn.c = c;
modelIn.X0 = X0;
%modelIn.InputParams.inputTransitionIdx = [1, 2];
modelIn.ObservedSpeciesIdx = 2;

if (generate == 1)
    [odeHandle, infos] = GenerateMomentsInputLN(modelIn);
    save odeInfosLN.mat odeHandle infos;
else
    load odeInfosLN.mat;
end


