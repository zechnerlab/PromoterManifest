function [dM, Likelihood] = ODEInputMomentsLN(t, M, model)

dM = M;
update=model.Update;
c2 = model.c(2);
c3 = model.c(3);
c4 = model.c(4);
X_1 = M(1);
X_2 = M(2);
X_3 = M(3);
X_11 = M(4);
X_12 = M(5);
X_13 = M(6);
X_22 = M(7);
X_23 = M(8);
X_33 = M(9);
X_113 = M(10);
X_123 = M(11);
X_133 = M(12);
X_1133 = M(13);
Z = model.Z;
if (update==0)
	dM(1) = Z - X_1*c2;
	dM(2) = X_13*c3 - X_2*c4;
	dM(3) = 0;
	dM(4) = Z + 2*X_1*Z + X_1*c2 - 2*X_11*c2;
	dM(5) = X_2*Z - X_12*c2 - X_12*c4 + X_113*c3;
	dM(6) = X_3*Z - X_13*c2;
	dM(7) = X_2*c4 + X_13*c3 - 2*X_22*c4 + 2*X_123*c3;
	dM(8) = X_133*c3 - X_23*c4;
	dM(9) = 0;
	dM(10) = X_3*Z + 2*X_13*Z + X_13*c2 - 2*X_113*c2;
	dM(11) = X_23*Z - X_123*c2 - X_123*c4 + X_1133*c3;
	dM(12) = X_33*Z - X_133*c2;
	dM(13) = X_33*Z + 2*X_133*Z + X_133*c2 - 2*X_1133*c2;
else
    warning off;
	w = zeros(3, 1);
	w(2) = 1;
	E=ones(3, 1);
	v=model.MeasurementSigma;
    
	y=log(model.Measurement+10); %add a small numerical value to avoid 0 inside the log. 

	MuPrior = zeros(3, 1);
	MuPrior(1) = X_1;
	MuPrior(2) = X_2;
	MuPrior(3) = X_3;
	SPrior(1, 1) = X_11;
	SPrior(1, 1) = X_11;
	SPrior(1, 2) = X_12;
	SPrior(2, 1) = X_12;
	SPrior(1, 3) = X_13;
	SPrior(3, 1) = X_13;
	SPrior(2, 2) = X_22;
	SPrior(2, 2) = X_22;
	SPrior(2, 3) = X_23;
	SPrior(3, 2) = X_23;
	SPrior(3, 3) = X_33;
	SPrior(3, 3) = X_33;
	KPrior(1) = X_113;
	KPrior(2) = X_123;
	KPrior(3) = X_133;
	UPrior(1) = X_1133;
	MuPriorLN = 2*log(MuPrior) - 1/2*log(diag(SPrior));
    %SigmaPriorLN = log(SPrior);
    
    
    
	SigmaPriorLN = log(SPrior) - log(MuPrior)*E' - E*log(MuPrior)';
	InvSigmaPriorLN = inv(SigmaPriorLN);
	InvSigmaLN = w*w' / v^2 + InvSigmaPriorLN;
	SigmaPosteriorLN = inv(InvSigmaLN);
	%MuPosteriorLN = invSigmaLN\(w*y / v^2 + SigmaPriorLN\MuPriorLN);
    %MuPosteriorLN = InvSigmaLN\(w*y / v^2 + (SigmaPriorLN\MuPriorLN));
    MuPosteriorLN = SigmaPosteriorLN * (w*y / v^2 + InvSigmaPriorLN * MuPriorLN);
    
	m = w'*MuPriorLN;
	s = SigmaPriorLN(2, 2);
	%Likelihood = -((m^2 + y^2)/(2*(s + v^2))) + (-1 + m/(s + v^2))*y - 1/2*log(s + v^2);
    Likelihood = -1/2*(y-m)^2/(v^2 + s) - 1/2*log(1/v^2 + 1/s) - 1/2*log(s) - log(v);
    
%     if sum(diag(SigmaPosteriorLN)<0)>0
%        fprintf('nan\n'); 
%     end
%     if (abs(imag(Likelihood))>0)
%        fprintf('nan\n'); 
%     end
% 	if isnan(Likelihood) || isinf(Likelihood) || abs(imag(Likelihood))>0
% 		Likelihood = -inf;
% 		MuPosterior = zeros(size(MuPosterior));
% 		SPosterior = zeros(size(SPosterior));
% 		SigmaPosterior = zeros(size(SigmaPosterior));
% 	else
% 	end
	MuPosteriorZ=exp(MuPosteriorLN + 1/2*diag(SigmaPosteriorLN));
	SPosteriorZ(1, 1) = exp(MuPosteriorLN(1) + MuPosteriorLN(1) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 1)+...
		SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 1)));
	SPosteriorZ(1, 1) = SPosteriorZ(1, 1);
	SPosteriorZ(1, 2) = exp(MuPosteriorLN(1) + MuPosteriorLN(2) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 2)+...
		SigmaPosteriorLN(2, 1) + SigmaPosteriorLN(2, 2)));
	SPosteriorZ(2, 1) = SPosteriorZ(1, 2);
	SPosteriorZ(1, 3) = exp(MuPosteriorLN(1) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3)+...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3)));
	SPosteriorZ(3, 1) = SPosteriorZ(1, 3);
	SPosteriorZ(2, 2) = exp(MuPosteriorLN(2) + MuPosteriorLN(2) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 2)+...
		SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 2)));
	SPosteriorZ(2, 2) = SPosteriorZ(2, 2);
	SPosteriorZ(2, 3) = exp(MuPosteriorLN(2) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 3)+...
		SigmaPosteriorLN(3, 2) + SigmaPosteriorLN(3, 3)));
	SPosteriorZ(3, 2) = SPosteriorZ(2, 3);
	SPosteriorZ(3, 3) = exp(MuPosteriorLN(3) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3)+...
		SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3)));
	SPosteriorZ(3, 3) = SPosteriorZ(3, 3);

	MuPosterior=MuPosteriorZ;
	SPosterior=SPosteriorZ;

	KPosterior(1) = exp(MuPosteriorLN(1) + MuPosteriorLN(1) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3) +...
		SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3) + ...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3)));

	KPosterior(2) = exp(MuPosteriorLN(1) + MuPosteriorLN(2) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 2) + SigmaPosteriorLN(1, 3) +...
		SigmaPosteriorLN(2, 1) + SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 3) + ...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 2) + SigmaPosteriorLN(3, 3)));

	KPosterior(3) = exp(MuPosteriorLN(1) + MuPosteriorLN(3) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3) + SigmaPosteriorLN(1, 3) +...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3) + ...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3)));

	UPosterior(1) = exp(MuPosteriorLN(1) + MuPosteriorLN(1) + MuPosteriorLN(3) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3) + SigmaPosteriorLN(1, 3) + ...
		SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3) + SigmaPosteriorLN(1, 3) +...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3) +...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3)));

	% Store values in change vector
	dM(1) = MuPosterior(1) - MuPrior(1);
	dM(2) = MuPosterior(2) - MuPrior(2);
	dM(3) = MuPosterior(3) - MuPrior(3);
	dM(4) = SPosterior(1, 1) - SPrior(1, 1);
	dM(5) = SPosterior(1, 2) - SPrior(1, 2);
	dM(6) = SPosterior(1, 3) - SPrior(1, 3);
	dM(7) = SPosterior(2, 2) - SPrior(2, 2);
	dM(8) = SPosterior(2, 3) - SPrior(2, 3);
	dM(9) = SPosterior(3, 3) - SPrior(3, 3);
	dM(10) = KPosterior(1)- KPrior(1);
	dM(11) = KPosterior(2)- KPrior(2);
	dM(12) = KPosterior(3)- KPrior(3);
	dM(13) = UPosterior(1) - UPrior(1);
end
