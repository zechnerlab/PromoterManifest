clear all;
close all;

addpath('../Common');



[Pre, Post, c, X0] = CreateTwoStageExpressionSystemZ();
X0(1) = 0;
X0(2) = 0;
%X0(3) = 0;


modelIn.Pre = Pre;
modelIn.Post = Post;
modelIn.c = c;
modelIn.X0 = X0;
modelIn.InputParams.inputTransitionIdx = [1, 2];
modelIn.ObservedSpeciesIdx = 2;

numStates = 3;


[odeHandle, infos] = GenerateConditionalBackwardMomentsLN(modelIn, numStates);
save odeInfosBackwardLN.mat odeHandle infos;
