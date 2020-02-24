function [dM, Likelihood, LTot] = ODECondMoments(t, M, model)

dM = M;
update=model.Update;
c2 = model.c(2);
c3 = model.c(3);
c4 = model.c(4);
z_1 = model.Z(1);
z_2 = model.Z(2);
z_3 = model.Z(3);
P = M(1:3);
X_1_1 = M(4);
X_1_2 = M(5);
X_1_3 = M(6);
X_1_11 = M(7);
X_1_12 = M(8);
X_1_13 = M(9);
X_1_22 = M(10);
X_1_23 = M(11);
X_1_33 = M(12);
X_1_113 = M(13);
X_1_123 = M(14);
X_1_133 = M(15);
X_1_1133 = M(16);
X_2_1 = M(17);
X_2_2 = M(18);
X_2_3 = M(19);
X_2_11 = M(20);
X_2_12 = M(21);
X_2_13 = M(22);
X_2_22 = M(23);
X_2_23 = M(24);
X_2_33 = M(25);
X_2_113 = M(26);
X_2_123 = M(27);
X_2_133 = M(28);
X_2_1133 = M(29);
X_3_1 = M(30);
X_3_2 = M(31);
X_3_3 = M(32);
X_3_11 = M(33);
X_3_12 = M(34);
X_3_13 = M(35);
X_3_22 = M(36);
X_3_23 = M(37);
X_3_33 = M(38);
X_3_113 = M(39);
X_3_123 = M(40);
X_3_133 = M(41);
X_3_1133 = M(42);
u = GetPieceWiseConstantInput(t, model.InputParams);q_2_1 = model.Q(1, 2);
q_3_1 = model.Q(1, 3);
q_1_2 = model.Q(2, 1);
q_3_2 = model.Q(2, 3);
q_1_3 = model.Q(3, 1);
q_2_3 = model.Q(3, 2);

if (update==0)
	dM(1) = -P(1)*(q_1_3 + q_1_2*u)+ q_2_1*P(2)+ q_3_1*P(3);
	dM(2) = q_1_2*u*P(1)+ -P(2)*(q_2_1 + q_2_3)+ q_3_2*P(3);
	dM(3) = q_1_3*P(1)+ q_2_3*P(2)+ -P(3)*(q_3_1 + q_3_2);
	dM(4) = X_2_1*q_2_1 - X_1_1*q_1_3 - X_1_1*c2 + X_3_1*q_3_1 + z_1*P(1) - X_1_1*q_1_2*u;
	dM(5) = X_1_13*c3 - X_1_2*c4 - X_1_2*q_1_3 + X_2_2*q_2_1 + X_3_2*q_3_1 - X_1_2*q_1_2*u;
	dM(6) = X_2_3*q_2_1 - X_1_3*q_1_3 + X_3_3*q_3_1 - X_1_3*q_1_2*u;
	dM(7) = z_1 + X_1_1*c2 - 2*X_1_11*c2 - X_1_11*q_1_3 + X_2_11*q_2_1 + X_3_11*q_3_1 + 2*X_1_1*z_1 - X_1_11*q_1_2*u;
	dM(8) = X_1_113*c3 - X_1_12*c4 - X_1_12*c2 - X_1_12*q_1_3 + X_2_12*q_2_1 + X_3_12*q_3_1 + X_1_2*z_1 - X_1_12*q_1_2*u;
	dM(9) = X_2_13*q_2_1 - X_1_13*q_1_3 - X_1_13*c2 + X_3_13*q_3_1 + X_1_3*z_1 - X_1_13*q_1_2*u;
	dM(10) = X_1_2*c4 + X_1_13*c3 - 2*X_1_22*c4 + 2*X_1_123*c3 - X_1_22*q_1_3 + X_2_22*q_2_1 + X_3_22*q_3_1 - X_1_22*q_1_2*u;
	dM(11) = X_1_133*c3 - X_1_23*c4 - X_1_23*q_1_3 + X_2_23*q_2_1 + X_3_23*q_3_1 - X_1_23*q_1_2*u;
	dM(12) = X_2_33*q_2_1 - X_1_33*q_1_3 + X_3_33*q_3_1 - X_1_33*q_1_2*u;
	dM(13) = X_1_13*c2 - 2*X_1_113*c2 - X_1_113*q_1_3 + X_2_113*q_2_1 + X_3_113*q_3_1 + X_1_3*z_1 + 2*X_1_13*z_1 - X_1_113*q_1_2*u;
	dM(14) = X_1_1133*c3 - X_1_123*c4 - X_1_123*c2 - X_1_123*q_1_3 + X_2_123*q_2_1 + X_3_123*q_3_1 + X_1_23*z_1 - X_1_123*q_1_2*u;
	dM(15) = X_2_133*q_2_1 - X_1_133*q_1_3 - X_1_133*c2 + X_3_133*q_3_1 + X_1_33*z_1 - X_1_133*q_1_2*u;
	dM(16) = X_1_133*c2 - 2*X_1_1133*c2 - X_1_1133*q_1_3 + X_2_1133*q_2_1 + X_3_1133*q_3_1 + X_1_33*z_1 + 2*X_1_133*z_1 - X_1_1133*q_1_2*u;
	dM(17) = X_3_1*q_3_2 - X_2_1*q_2_1 - X_2_1*q_2_3 - X_2_1*c2 + z_2*P(2) + X_1_1*q_1_2*u;
	dM(18) = X_2_13*c3 - X_2_2*c4 - X_2_2*q_2_1 - X_2_2*q_2_3 + X_3_2*q_3_2 + X_1_2*q_1_2*u;
	dM(19) = X_3_3*q_3_2 - X_2_3*q_2_3 - X_2_3*q_2_1 + X_1_3*q_1_2*u;
	dM(20) = z_2 + X_2_1*c2 - 2*X_2_11*c2 - X_2_11*q_2_1 - X_2_11*q_2_3 + X_3_11*q_3_2 + 2*X_2_1*z_2 + X_1_11*q_1_2*u;
	dM(21) = X_2_113*c3 - X_2_12*c4 - X_2_12*c2 - X_2_12*q_2_1 - X_2_12*q_2_3 + X_3_12*q_3_2 + X_2_2*z_2 + X_1_12*q_1_2*u;
	dM(22) = X_3_13*q_3_2 - X_2_13*q_2_1 - X_2_13*q_2_3 - X_2_13*c2 + X_2_3*z_2 + X_1_13*q_1_2*u;
	dM(23) = X_2_2*c4 + X_2_13*c3 - 2*X_2_22*c4 + 2*X_2_123*c3 - X_2_22*q_2_1 - X_2_22*q_2_3 + X_3_22*q_3_2 + X_1_22*q_1_2*u;
	dM(24) = X_2_133*c3 - X_2_23*c4 - X_2_23*q_2_1 - X_2_23*q_2_3 + X_3_23*q_3_2 + X_1_23*q_1_2*u;
	dM(25) = X_3_33*q_3_2 - X_2_33*q_2_3 - X_2_33*q_2_1 + X_1_33*q_1_2*u;
	dM(26) = X_2_13*c2 - 2*X_2_113*c2 - X_2_113*q_2_1 - X_2_113*q_2_3 + X_3_113*q_3_2 + X_2_3*z_2 + 2*X_2_13*z_2 + X_1_113*q_1_2*u;
	dM(27) = X_2_1133*c3 - X_2_123*c4 - X_2_123*c2 - X_2_123*q_2_1 - X_2_123*q_2_3 + X_3_123*q_3_2 + X_2_23*z_2 + X_1_123*q_1_2*u;
	dM(28) = X_3_133*q_3_2 - X_2_133*q_2_1 - X_2_133*q_2_3 - X_2_133*c2 + X_2_33*z_2 + X_1_133*q_1_2*u;
	dM(29) = X_2_133*c2 - 2*X_2_1133*c2 - X_2_1133*q_2_1 - X_2_1133*q_2_3 + X_3_1133*q_3_2 + X_2_33*z_2 + 2*X_2_133*z_2 + X_1_1133*q_1_2*u;
	dM(30) = X_1_1*q_1_3 - X_3_1*c2 + X_2_1*q_2_3 - X_3_1*q_3_1 - X_3_1*q_3_2 + z_3*P(3);
	dM(31) = X_3_13*c3 - X_3_2*c4 + X_1_2*q_1_3 + X_2_2*q_2_3 - X_3_2*q_3_1 - X_3_2*q_3_2;
	dM(32) = X_1_3*q_1_3 + X_2_3*q_2_3 - X_3_3*q_3_1 - X_3_3*q_3_2;
	dM(33) = z_3 + X_3_1*c2 - 2*X_3_11*c2 + X_1_11*q_1_3 + X_2_11*q_2_3 - X_3_11*q_3_1 - X_3_11*q_3_2 + 2*X_3_1*z_3;
	dM(34) = X_3_113*c3 - X_3_12*c4 - X_3_12*c2 + X_1_12*q_1_3 + X_2_12*q_2_3 - X_3_12*q_3_1 - X_3_12*q_3_2 + X_3_2*z_3;
	dM(35) = X_1_13*q_1_3 - X_3_13*c2 + X_2_13*q_2_3 - X_3_13*q_3_1 - X_3_13*q_3_2 + X_3_3*z_3;
	dM(36) = X_3_2*c4 + X_3_13*c3 - 2*X_3_22*c4 + 2*X_3_123*c3 + X_1_22*q_1_3 + X_2_22*q_2_3 - X_3_22*q_3_1 - X_3_22*q_3_2;
	dM(37) = X_3_133*c3 - X_3_23*c4 + X_1_23*q_1_3 + X_2_23*q_2_3 - X_3_23*q_3_1 - X_3_23*q_3_2;
	dM(38) = X_1_33*q_1_3 + X_2_33*q_2_3 - X_3_33*q_3_1 - X_3_33*q_3_2;
	dM(39) = X_3_13*c2 - 2*X_3_113*c2 + X_1_113*q_1_3 + X_2_113*q_2_3 - X_3_113*q_3_1 - X_3_113*q_3_2 + X_3_3*z_3 + 2*X_3_13*z_3;
	dM(40) = X_3_1133*c3 - X_3_123*c4 - X_3_123*c2 + X_1_123*q_1_3 + X_2_123*q_2_3 - X_3_123*q_3_1 - X_3_123*q_3_2 + X_3_23*z_3;
	dM(41) = X_1_133*q_1_3 - X_3_133*c2 + X_2_133*q_2_3 - X_3_133*q_3_1 - X_3_133*q_3_2 + X_3_33*z_3;
	dM(42) = X_3_133*c2 - 2*X_3_1133*c2 + X_1_1133*q_1_3 + X_2_1133*q_2_3 - X_3_1133*q_3_1 - X_3_1133*q_3_2 + X_3_33*z_3 + 2*X_3_133*z_3;
else
	PPost = zeros(1, 3);
	w = zeros(3, 1);
	w(2) = 1;
	E=ones(3, 1);
	v=model.MeasurementSigma;
	y=log(model.Measurement+1);

	%State 1 
	MuPrior1 = zeros(3, 1);
	MuPrior1(1) = X_1_1 + eps;
	MuPrior1(2) = X_1_2 + eps;
	MuPrior1(3) = X_1_3 + eps;
	SPrior1(1, 1) = X_1_11 + eps;
	SPrior1(1, 1) = X_1_11 + eps;
	SPrior1(1, 2) = X_1_12 + eps;
	SPrior1(2, 1) = X_1_12 + eps;
	SPrior1(1, 3) = X_1_13 + eps;
	SPrior1(3, 1) = X_1_13 + eps;
	SPrior1(2, 2) = X_1_22 + eps;
	SPrior1(2, 2) = X_1_22 + eps;
	SPrior1(2, 3) = X_1_23 + eps;
	SPrior1(3, 2) = X_1_23 + eps;
	SPrior1(3, 3) = X_1_33 + eps;
	SPrior1(3, 3) = X_1_33 + eps;
	KPrior1(1) = X_1_113 + eps;
	KPrior1(2) = X_1_123 + eps;
	KPrior1(3) = X_1_133 + eps;
	UPrior1(1) = X_1_1133 + eps;
	MuPriorLN1 = 2*log(MuPrior1) - 3/2*log(P(1)) - 1/2*log(diag(SPrior1));
	SigmaPriorLN1 = -log(P(1)) + log(SPrior1) - log(MuPrior1)*E' - E*log(MuPrior1)' + 2*log(P(1));
	InvSigmaPriorLN = inv(SigmaPriorLN1);
	InvSigmaLN = w*w' / v^2 + InvSigmaPriorLN;
	SigmaPosteriorLN1 = inv(InvSigmaLN);
	MuPosteriorLN1 = SigmaPosteriorLN1 * (w*y / v^2 + InvSigmaPriorLN * MuPriorLN1);
	m = w'*MuPriorLN1;
	s = SigmaPriorLN1(2, 2);
	Likelihood(1) = -((m^2 + y^2)/(2*(s + v^2))) + (-1 + m/(s + v^2))*y - 1/2*log(s + v^2);
	if isnan(Likelihood(1)) || isinf(Likelihood(1))
		Likelihood(1) = -inf;
		PPost(1) = 0;
		MuPosterior1 = zeros(size(MuPosterior1));
		SPosterior1 = zeros(size(SPosterior1));
		SigmaPosterior1 = zeros(size(SigmaPosterior1));
	else
		PPost(1) = exp(Likelihood(1))*P(1);
	end
	MuPosteriorZ1=exp(MuPosteriorLN1 + 1/2*diag(SigmaPosteriorLN1));
	SPosteriorZ1(1, 1) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(1) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 1)+...
		SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 1)));
	SPosteriorZ1(1, 1) = SPosteriorZ1(1, 1);
	SPosteriorZ1(1, 2) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(2) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 2)+...
		SigmaPosteriorLN1(2, 1) + SigmaPosteriorLN1(2, 2)));
	SPosteriorZ1(2, 1) = SPosteriorZ1(1, 2);
	SPosteriorZ1(1, 3) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 3)+...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 3)));
	SPosteriorZ1(3, 1) = SPosteriorZ1(1, 3);
	SPosteriorZ1(2, 2) = exp(MuPosteriorLN1(2) + MuPosteriorLN1(2) + ...
		1/2*(SigmaPosteriorLN1(2, 2) + SigmaPosteriorLN1(2, 2)+...
		SigmaPosteriorLN1(2, 2) + SigmaPosteriorLN1(2, 2)));
	SPosteriorZ1(2, 2) = SPosteriorZ1(2, 2);
	SPosteriorZ1(2, 3) = exp(MuPosteriorLN1(2) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(2, 2) + SigmaPosteriorLN1(2, 3)+...
		SigmaPosteriorLN1(3, 2) + SigmaPosteriorLN1(3, 3)));
	SPosteriorZ1(3, 2) = SPosteriorZ1(2, 3);
	SPosteriorZ1(3, 3) = exp(MuPosteriorLN1(3) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(3, 3) + SigmaPosteriorLN1(3, 3)+...
		SigmaPosteriorLN1(3, 3) + SigmaPosteriorLN1(3, 3)));
	SPosteriorZ1(3, 3) = SPosteriorZ1(3, 3);

	MuPosterior1=MuPosteriorZ1;
	SPosterior1=SPosteriorZ1;

	KPosterior1(1) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(1) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 3) +...
		SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 3) + ...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 3)));

	KPosterior1(2) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(2) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 2) + SigmaPosteriorLN1(1, 3) +...
		SigmaPosteriorLN1(2, 1) + SigmaPosteriorLN1(2, 2) + SigmaPosteriorLN1(2, 3) + ...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 2) + SigmaPosteriorLN1(3, 3)));

	KPosterior1(3) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(3) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 3) + SigmaPosteriorLN1(1, 3) +...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 3) + SigmaPosteriorLN1(3, 3) + ...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 3) + SigmaPosteriorLN1(3, 3)));

	UPosterior1(1) = exp(MuPosteriorLN1(1) + MuPosteriorLN1(1) + MuPosteriorLN1(3) + MuPosteriorLN1(3) + ...
		1/2*(SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 3) + SigmaPosteriorLN1(1, 3) + ...
		SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 1) + SigmaPosteriorLN1(1, 3) + SigmaPosteriorLN1(1, 3) +...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 3) + SigmaPosteriorLN1(3, 3) +...
		SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 1) + SigmaPosteriorLN1(3, 3) + SigmaPosteriorLN1(3, 3)));


	%State 2 
	MuPrior2 = zeros(3, 1);
	MuPrior2(1) = X_2_1 + eps;
	MuPrior2(2) = X_2_2 + eps;
	MuPrior2(3) = X_2_3 + eps;
	SPrior2(1, 1) = X_2_11 + eps;
	SPrior2(1, 1) = X_2_11 + eps;
	SPrior2(1, 2) = X_2_12 + eps;
	SPrior2(2, 1) = X_2_12 + eps;
	SPrior2(1, 3) = X_2_13 + eps;
	SPrior2(3, 1) = X_2_13 + eps;
	SPrior2(2, 2) = X_2_22 + eps;
	SPrior2(2, 2) = X_2_22 + eps;
	SPrior2(2, 3) = X_2_23 + eps;
	SPrior2(3, 2) = X_2_23 + eps;
	SPrior2(3, 3) = X_2_33 + eps;
	SPrior2(3, 3) = X_2_33 + eps;
	KPrior2(1) = X_2_113 + eps;
	KPrior2(2) = X_2_123 + eps;
	KPrior2(3) = X_2_133 + eps;
	UPrior2(1) = X_2_1133 + eps;
	MuPriorLN2 = 2*log(MuPrior2) - 3/2*log(P(2)) - 1/2*log(diag(SPrior2));
	SigmaPriorLN2 = -log(P(2)) + log(SPrior2) - log(MuPrior2)*E' - E*log(MuPrior2)' + 2*log(P(2));
	InvSigmaPriorLN = inv(SigmaPriorLN2);
	InvSigmaLN = w*w' / v^2 + InvSigmaPriorLN;
	SigmaPosteriorLN2 = inv(InvSigmaLN);
	MuPosteriorLN2 = SigmaPosteriorLN2 * (w*y / v^2 + InvSigmaPriorLN * MuPriorLN2);
	m = w'*MuPriorLN2;
	s = SigmaPriorLN2(2, 2);
	Likelihood(2) = -((m^2 + y^2)/(2*(s + v^2))) + (-1 + m/(s + v^2))*y - 1/2*log(s + v^2);
	if isnan(Likelihood(2)) || isinf(Likelihood(2))
		Likelihood(2) = -inf;
		PPost(2) = 0;
		MuPosterior2 = zeros(size(MuPosterior2));
		SPosterior2 = zeros(size(SPosterior2));
		SigmaPosterior2 = zeros(size(SigmaPosterior2));
	else
		PPost(2) = exp(Likelihood(2))*P(2);
	end
	MuPosteriorZ2=exp(MuPosteriorLN2 + 1/2*diag(SigmaPosteriorLN2));
	SPosteriorZ2(1, 1) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(1) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 1)+...
		SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 1)));
	SPosteriorZ2(1, 1) = SPosteriorZ2(1, 1);
	SPosteriorZ2(1, 2) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(2) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 2)+...
		SigmaPosteriorLN2(2, 1) + SigmaPosteriorLN2(2, 2)));
	SPosteriorZ2(2, 1) = SPosteriorZ2(1, 2);
	SPosteriorZ2(1, 3) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 3)+...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 3)));
	SPosteriorZ2(3, 1) = SPosteriorZ2(1, 3);
	SPosteriorZ2(2, 2) = exp(MuPosteriorLN2(2) + MuPosteriorLN2(2) + ...
		1/2*(SigmaPosteriorLN2(2, 2) + SigmaPosteriorLN2(2, 2)+...
		SigmaPosteriorLN2(2, 2) + SigmaPosteriorLN2(2, 2)));
	SPosteriorZ2(2, 2) = SPosteriorZ2(2, 2);
	SPosteriorZ2(2, 3) = exp(MuPosteriorLN2(2) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(2, 2) + SigmaPosteriorLN2(2, 3)+...
		SigmaPosteriorLN2(3, 2) + SigmaPosteriorLN2(3, 3)));
	SPosteriorZ2(3, 2) = SPosteriorZ2(2, 3);
	SPosteriorZ2(3, 3) = exp(MuPosteriorLN2(3) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(3, 3) + SigmaPosteriorLN2(3, 3)+...
		SigmaPosteriorLN2(3, 3) + SigmaPosteriorLN2(3, 3)));
	SPosteriorZ2(3, 3) = SPosteriorZ2(3, 3);

	MuPosterior2=MuPosteriorZ2;
	SPosterior2=SPosteriorZ2;

	KPosterior2(1) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(1) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 3) +...
		SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 3) + ...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 3)));

	KPosterior2(2) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(2) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 2) + SigmaPosteriorLN2(1, 3) +...
		SigmaPosteriorLN2(2, 1) + SigmaPosteriorLN2(2, 2) + SigmaPosteriorLN2(2, 3) + ...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 2) + SigmaPosteriorLN2(3, 3)));

	KPosterior2(3) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(3) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 3) + SigmaPosteriorLN2(1, 3) +...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 3) + SigmaPosteriorLN2(3, 3) + ...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 3) + SigmaPosteriorLN2(3, 3)));

	UPosterior2(1) = exp(MuPosteriorLN2(1) + MuPosteriorLN2(1) + MuPosteriorLN2(3) + MuPosteriorLN2(3) + ...
		1/2*(SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 3) + SigmaPosteriorLN2(1, 3) + ...
		SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 1) + SigmaPosteriorLN2(1, 3) + SigmaPosteriorLN2(1, 3) +...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 3) + SigmaPosteriorLN2(3, 3) +...
		SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 1) + SigmaPosteriorLN2(3, 3) + SigmaPosteriorLN2(3, 3)));


	%State 3 
	MuPrior3 = zeros(3, 1);
	MuPrior3(1) = X_3_1 + eps;
	MuPrior3(2) = X_3_2 + eps;
	MuPrior3(3) = X_3_3 + eps;
	SPrior3(1, 1) = X_3_11 + eps;
	SPrior3(1, 1) = X_3_11 + eps;
	SPrior3(1, 2) = X_3_12 + eps;
	SPrior3(2, 1) = X_3_12 + eps;
	SPrior3(1, 3) = X_3_13 + eps;
	SPrior3(3, 1) = X_3_13 + eps;
	SPrior3(2, 2) = X_3_22 + eps;
	SPrior3(2, 2) = X_3_22 + eps;
	SPrior3(2, 3) = X_3_23 + eps;
	SPrior3(3, 2) = X_3_23 + eps;
	SPrior3(3, 3) = X_3_33 + eps;
	SPrior3(3, 3) = X_3_33 + eps;
	KPrior3(1) = X_3_113 + eps;
	KPrior3(2) = X_3_123 + eps;
	KPrior3(3) = X_3_133 + eps;
	UPrior3(1) = X_3_1133 + eps;
	MuPriorLN3 = 2*log(MuPrior3) - 3/2*log(P(3)) - 1/2*log(diag(SPrior3));
	SigmaPriorLN3 = -log(P(3)) + log(SPrior3) - log(MuPrior3)*E' - E*log(MuPrior3)' + 2*log(P(3));
	InvSigmaPriorLN = inv(SigmaPriorLN3);
	InvSigmaLN = w*w' / v^2 + InvSigmaPriorLN;
	SigmaPosteriorLN3 = inv(InvSigmaLN);
	MuPosteriorLN3 = SigmaPosteriorLN3 * (w*y / v^2 + InvSigmaPriorLN * MuPriorLN3);
	m = w'*MuPriorLN3;
	s = SigmaPriorLN3(2, 2);
	Likelihood(3) = -((m^2 + y^2)/(2*(s + v^2))) + (-1 + m/(s + v^2))*y - 1/2*log(s + v^2);
	if isnan(Likelihood(3)) || isinf(Likelihood(3))
		Likelihood(3) = -inf;
		PPost(3) = 0;
		MuPosterior3 = zeros(size(MuPosterior3));
		SPosterior3 = zeros(size(SPosterior3));
		SigmaPosterior3 = zeros(size(SigmaPosterior3));
	else
		PPost(3) = exp(Likelihood(3))*P(3);
	end
	MuPosteriorZ3=exp(MuPosteriorLN3 + 1/2*diag(SigmaPosteriorLN3));
	SPosteriorZ3(1, 1) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(1) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 1)+...
		SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 1)));
	SPosteriorZ3(1, 1) = SPosteriorZ3(1, 1);
	SPosteriorZ3(1, 2) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(2) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 2)+...
		SigmaPosteriorLN3(2, 1) + SigmaPosteriorLN3(2, 2)));
	SPosteriorZ3(2, 1) = SPosteriorZ3(1, 2);
	SPosteriorZ3(1, 3) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 3)+...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 3)));
	SPosteriorZ3(3, 1) = SPosteriorZ3(1, 3);
	SPosteriorZ3(2, 2) = exp(MuPosteriorLN3(2) + MuPosteriorLN3(2) + ...
		1/2*(SigmaPosteriorLN3(2, 2) + SigmaPosteriorLN3(2, 2)+...
		SigmaPosteriorLN3(2, 2) + SigmaPosteriorLN3(2, 2)));
	SPosteriorZ3(2, 2) = SPosteriorZ3(2, 2);
	SPosteriorZ3(2, 3) = exp(MuPosteriorLN3(2) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(2, 2) + SigmaPosteriorLN3(2, 3)+...
		SigmaPosteriorLN3(3, 2) + SigmaPosteriorLN3(3, 3)));
	SPosteriorZ3(3, 2) = SPosteriorZ3(2, 3);
	SPosteriorZ3(3, 3) = exp(MuPosteriorLN3(3) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(3, 3) + SigmaPosteriorLN3(3, 3)+...
		SigmaPosteriorLN3(3, 3) + SigmaPosteriorLN3(3, 3)));
	SPosteriorZ3(3, 3) = SPosteriorZ3(3, 3);

	MuPosterior3=MuPosteriorZ3;
	SPosterior3=SPosteriorZ3;

	KPosterior3(1) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(1) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 3) +...
		SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 3) + ...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 3)));

	KPosterior3(2) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(2) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 2) + SigmaPosteriorLN3(1, 3) +...
		SigmaPosteriorLN3(2, 1) + SigmaPosteriorLN3(2, 2) + SigmaPosteriorLN3(2, 3) + ...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 2) + SigmaPosteriorLN3(3, 3)));

	KPosterior3(3) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(3) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 3) + SigmaPosteriorLN3(1, 3) +...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 3) + SigmaPosteriorLN3(3, 3) + ...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 3) + SigmaPosteriorLN3(3, 3)));

	UPosterior3(1) = exp(MuPosteriorLN3(1) + MuPosteriorLN3(1) + MuPosteriorLN3(3) + MuPosteriorLN3(3) + ...
		1/2*(SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 3) + SigmaPosteriorLN3(1, 3) + ...
		SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 1) + SigmaPosteriorLN3(1, 3) + SigmaPosteriorLN3(1, 3) +...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 3) + SigmaPosteriorLN3(3, 3) +...
		SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 1) + SigmaPosteriorLN3(3, 3) + SigmaPosteriorLN3(3, 3)));


	PPost = PPost / sum(PPost);

	% Store values in change vector
	dM(1:3) = PPost - P;
	dM(4) = MuPosterior1(1)*PPost(1) - MuPrior1(1);
	dM(5) = MuPosterior1(2)*PPost(1) - MuPrior1(2);
	dM(6) = MuPosterior1(3)*PPost(1) - MuPrior1(3);
	dM(7) = SPosterior1(1, 1)*PPost(1) - SPrior1(1, 1);
	dM(8) = SPosterior1(1, 2)*PPost(1) - SPrior1(1, 2);
	dM(9) = SPosterior1(1, 3)*PPost(1) - SPrior1(1, 3);
	dM(10) = SPosterior1(2, 2)*PPost(1) - SPrior1(2, 2);
	dM(11) = SPosterior1(2, 3)*PPost(1) - SPrior1(2, 3);
	dM(12) = SPosterior1(3, 3)*PPost(1) - SPrior1(3, 3);
	dM(13) = KPosterior1(1)*PPost(1) - KPrior1(1);
	dM(14) = KPosterior1(2)*PPost(1) - KPrior1(2);
	dM(15) = KPosterior1(3)*PPost(1) - KPrior1(3);
	dM(16) = UPosterior1(1)*PPost(1) - UPrior1(1);
	dM(17) = MuPosterior2(1)*PPost(2) - MuPrior2(1);
	dM(18) = MuPosterior2(2)*PPost(2) - MuPrior2(2);
	dM(19) = MuPosterior2(3)*PPost(2) - MuPrior2(3);
	dM(20) = SPosterior2(1, 1)*PPost(2) - SPrior2(1, 1);
	dM(21) = SPosterior2(1, 2)*PPost(2) - SPrior2(1, 2);
	dM(22) = SPosterior2(1, 3)*PPost(2) - SPrior2(1, 3);
	dM(23) = SPosterior2(2, 2)*PPost(2) - SPrior2(2, 2);
	dM(24) = SPosterior2(2, 3)*PPost(2) - SPrior2(2, 3);
	dM(25) = SPosterior2(3, 3)*PPost(2) - SPrior2(3, 3);
	dM(26) = KPosterior2(1)*PPost(2) - KPrior2(1);
	dM(27) = KPosterior2(2)*PPost(2) - KPrior2(2);
	dM(28) = KPosterior2(3)*PPost(2) - KPrior2(3);
	dM(29) = UPosterior2(1)*PPost(2) - UPrior2(1);
	dM(30) = MuPosterior3(1)*PPost(3) - MuPrior3(1);
	dM(31) = MuPosterior3(2)*PPost(3) - MuPrior3(2);
	dM(32) = MuPosterior3(3)*PPost(3) - MuPrior3(3);
	dM(33) = SPosterior3(1, 1)*PPost(3) - SPrior3(1, 1);
	dM(34) = SPosterior3(1, 2)*PPost(3) - SPrior3(1, 2);
	dM(35) = SPosterior3(1, 3)*PPost(3) - SPrior3(1, 3);
	dM(36) = SPosterior3(2, 2)*PPost(3) - SPrior3(2, 2);
	dM(37) = SPosterior3(2, 3)*PPost(3) - SPrior3(2, 3);
	dM(38) = SPosterior3(3, 3)*PPost(3) - SPrior3(3, 3);
	dM(39) = KPosterior3(1)*PPost(3) - KPrior3(1);
	dM(40) = KPosterior3(2)*PPost(3) - KPrior3(2);
	dM(41) = KPosterior3(3)*PPost(3) - KPrior3(3);
	dM(42) = UPosterior3(1)*PPost(3) - UPrior3(1);
	LTot = log(exp(Likelihood)*P');
	dM(isnan(dM)) = 0;
end
