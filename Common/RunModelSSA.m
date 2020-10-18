function [P, R, trProfiles] = RunModelSSA(model, time, inputParams, M)


[Pre, Post, c, X0] = CreateMSN2PromoterSystem();

c(1) = model.Q(2, 1);
c(2) = model.Q(1, 2);
c(3) = model.Q(3, 2);
c(4) = model.Q(2, 3);
c(5:7) = model.Z;
c(8:10) = model.c(2:end);
c(end) = 0.001;

alpha = 1/model.H(2);
beta = alpha;


inputParams.inputLevels = inputParams.inputLevels*c(1);

for k=1:M
   
   cTmp = c;
   Z = gamrnd(alpha, 1/beta);
   cTmp(9) = cTmp(9)*Z;
    

   [x, t] = SimulateSSA_TV(X0, cTmp, Pre, Post, 500000, time(end), 1, inputParams, [0, inputParams.inputTimes, time(end)], 0);
   
   Xs = SampleCTMPPathGrid_mex(x, t, time);
   
   P(k, :) = Xs(5, :);
   R(k, :) = Xs(4, :);
   trProfiles{k}.t = time;
   trProfiles{k}.Value = model.Z*Xs(1:3, :);
end

end