function [dM, Likelihood] = ODEInputMomentsLN(t, M, model)

dM = M;
update=model.Update;
c2 = model.c(2);
c3 = model.c(3);
c4 = model.c(4);
c5 = model.c(5);
c6 = model.c(6);
c7 = model.c(7);
c8 = model.c(8);
c9 = model.c(9);
X_1 = M(1);
X_2 = M(2);
X_3 = M(3);
X_4 = M(4);
X_5 = M(5);
X_6 = M(6);
X_11 = M(7);
X_12 = M(8);
X_13 = M(9);
X_14 = M(10);
X_15 = M(11);
X_16 = M(12);
X_22 = M(13);
X_23 = M(14);
X_24 = M(15);
X_25 = M(16);
X_26 = M(17);
X_33 = M(18);
X_34 = M(19);
X_35 = M(20);
X_36 = M(21);
X_44 = M(22);
X_45 = M(23);
X_46 = M(24);
X_55 = M(25);
X_56 = M(26);
X_66 = M(27);
X_146 = M(28);
X_246 = M(29);
X_346 = M(30);
X_446 = M(31);
X_456 = M(32);
X_466 = M(33);
X_4466 = M(34);
X_126 = M(35);
X_136 = M(36);
X_226 = M(37);
X_236 = M(38);
X_256 = M(39);
X_266 = M(40);
X_336 = M(41);
X_356 = M(42);
X_366 = M(43);
X_2466 = M(44);
X_3466 = M(45);
X_116 = M(46);
X_156 = M(47);
X_166 = M(48);
X_1466 = M(49);
X_2266 = M(50);
X_2366 = M(51);
X_3366 = M(52);
X_1266 = M(53);
X_1366 = M(54);
X_1166 = M(55);
Z = GetPieceWiseConstantInput(t, model.TranscriptionRateSequence);
if (update==0)
	dM(1) = X_2*c2 - X_1*Z;
	dM(2) = X_1*Z - X_2*c2 - X_2*c3 + X_3*c4;
	dM(3) = X_2*c3 - X_3*c4;
	dM(4) = X_2*c5 + X_3*c6 - X_4*c7;
	dM(5) = X_46*c8 - X_5*c9;
	dM(6) = 0;
	dM(7) = X_1*Z - 2*X_11*Z + X_2*c2 + 2*X_12*c2;
	dM(8) = X_11*Z - X_1*Z - X_12*Z - X_2*c2 - X_12*c2 - X_12*c3 + X_13*c4 + X_22*c2;
	dM(9) = X_12*c3 - X_13*Z - X_13*c4 + X_23*c2;
	dM(10) = X_12*c5 - X_14*Z + X_13*c6 - X_14*c7 + X_24*c2;
	dM(11) = X_25*c2 - X_15*c9 - X_15*Z + X_146*c8;
	dM(12) = X_26*c2 - X_16*Z;
	dM(13) = X_1*Z + 2*X_12*Z + X_2*c2 + X_2*c3 + X_3*c4 - 2*X_22*c2 - 2*X_22*c3 + 2*X_23*c4;
	dM(14) = X_13*Z - X_2*c3 - X_3*c4 + X_22*c3 - X_23*c2 - X_23*c3 - X_23*c4 + X_33*c4;
	dM(15) = X_14*Z - X_24*c2 + X_22*c5 - X_24*c3 + X_23*c6 - X_24*c7 + X_34*c4;
	dM(16) = X_15*Z - X_25*c2 - X_25*c3 - X_25*c9 + X_35*c4 + X_246*c8;
	dM(17) = X_16*Z - X_26*c2 - X_26*c3 + X_36*c4;
	dM(18) = X_2*c3 + X_3*c4 + 2*X_23*c3 - 2*X_33*c4;
	dM(19) = X_24*c3 + X_23*c5 - X_34*c4 + X_33*c6 - X_34*c7;
	dM(20) = X_25*c3 - X_35*c4 - X_35*c9 + X_346*c8;
	dM(21) = X_26*c3 - X_36*c4;
	dM(22) = X_2*c5 + X_3*c6 + X_4*c7 + 2*X_24*c5 + 2*X_34*c6 - 2*X_44*c7;
	dM(23) = X_25*c5 + X_35*c6 - X_45*c7 - X_45*c9 + X_446*c8;
	dM(24) = X_26*c5 + X_36*c6 - X_46*c7;
	dM(25) = X_5*c9 + X_46*c8 - 2*X_55*c9 + 2*X_456*c8;
	dM(26) = X_466*c8 - X_56*c9;
	dM(27) = 0;
	dM(28) = X_126*c5 - X_146*Z + X_136*c6 - X_146*c7 + X_246*c2;
	dM(29) = X_146*Z + X_226*c5 + X_236*c6 - X_246*c2 - X_246*c3 - X_246*c7 + X_346*c4;
	dM(30) = X_236*c5 + X_246*c3 + X_336*c6 - X_346*c4 - X_346*c7;
	dM(31) = X_26*c5 + X_36*c6 + X_46*c7 + 2*X_246*c5 + 2*X_346*c6 - 2*X_446*c7;
	dM(32) = X_256*c5 + X_356*c6 - X_456*c7 - X_456*c9 + X_4466*c8;
	dM(33) = X_266*c5 + X_366*c6 - X_466*c7;
	dM(34) = X_266*c5 + X_366*c6 + X_466*c7 + 2*X_2466*c5 + 2*X_3466*c6 - 2*X_4466*c7;
	dM(35) = X_116*Z - X_16*Z - X_126*Z - X_26*c2 - X_126*c2 - X_126*c3 + X_136*c4 + X_226*c2;
	dM(36) = X_126*c3 - X_136*Z - X_136*c4 + X_236*c2;
	dM(37) = X_16*Z + 2*X_126*Z + X_26*c2 + X_26*c3 + X_36*c4 - 2*X_226*c2 - 2*X_226*c3 + 2*X_236*c4;
	dM(38) = X_136*Z - X_26*c3 - X_36*c4 + X_226*c3 - X_236*c2 - X_236*c3 - X_236*c4 + X_336*c4;
	dM(39) = X_156*Z - X_256*c2 - X_256*c3 - X_256*c9 + X_356*c4 + X_2466*c8;
	dM(40) = X_166*Z - X_266*c2 - X_266*c3 + X_366*c4;
	dM(41) = X_26*c3 + X_36*c4 + 2*X_236*c3 - 2*X_336*c4;
	dM(42) = X_256*c3 - X_356*c4 - X_356*c9 + X_3466*c8;
	dM(43) = X_266*c3 - X_366*c4;
	dM(44) = X_1466*Z + X_2266*c5 + X_2366*c6 - X_2466*c2 - X_2466*c3 - X_2466*c7 + X_3466*c4;
	dM(45) = X_2366*c5 + X_2466*c3 + X_3366*c6 - X_3466*c4 - X_3466*c7;
	dM(46) = X_16*Z - 2*X_116*Z + X_26*c2 + 2*X_126*c2;
	dM(47) = X_256*c2 - X_156*c9 - X_156*Z + X_1466*c8;
	dM(48) = X_266*c2 - X_166*Z;
	dM(49) = X_1266*c5 - X_1466*Z + X_1366*c6 - X_1466*c7 + X_2466*c2;
	dM(50) = X_166*Z + 2*X_1266*Z + X_266*c2 + X_266*c3 + X_366*c4 - 2*X_2266*c2 - 2*X_2266*c3 + 2*X_2366*c4;
	dM(51) = X_1366*Z - X_266*c3 - X_366*c4 + X_2266*c3 - X_2366*c2 - X_2366*c3 - X_2366*c4 + X_3366*c4;
	dM(52) = X_266*c3 + X_366*c4 + 2*X_2366*c3 - 2*X_3366*c4;
	dM(53) = X_1166*Z - X_166*Z - X_1266*Z - X_266*c2 - X_1266*c2 - X_1266*c3 + X_1366*c4 + X_2266*c2;
	dM(54) = X_1266*c3 - X_1366*Z - X_1366*c4 + X_2366*c2;
	dM(55) = X_166*Z - 2*X_1166*Z + X_266*c2 + 2*X_1266*c2;
else
	w = zeros(6, 1);
	w(4) = 1;
	E=ones(6, 1);
	v=model.MeasurementSigma;
	y=log(model.Measurement+1);

	MuPrior = zeros(6, 1);
	MuPrior(1) = X_1;
	MuPrior(2) = X_2;
	MuPrior(3) = X_3;
	MuPrior(4) = X_4;
	MuPrior(5) = X_5;
	MuPrior(6) = X_6;
	SPrior(1, 1) = X_11;
	SPrior(1, 1) = X_11;
	SPrior(1, 2) = X_12;
	SPrior(2, 1) = X_12;
	SPrior(1, 3) = X_13;
	SPrior(3, 1) = X_13;
	SPrior(1, 4) = X_14;
	SPrior(4, 1) = X_14;
	SPrior(1, 5) = X_15;
	SPrior(5, 1) = X_15;
	SPrior(1, 6) = X_16;
	SPrior(6, 1) = X_16;
	SPrior(2, 2) = X_22;
	SPrior(2, 2) = X_22;
	SPrior(2, 3) = X_23;
	SPrior(3, 2) = X_23;
	SPrior(2, 4) = X_24;
	SPrior(4, 2) = X_24;
	SPrior(2, 5) = X_25;
	SPrior(5, 2) = X_25;
	SPrior(2, 6) = X_26;
	SPrior(6, 2) = X_26;
	SPrior(3, 3) = X_33;
	SPrior(3, 3) = X_33;
	SPrior(3, 4) = X_34;
	SPrior(4, 3) = X_34;
	SPrior(3, 5) = X_35;
	SPrior(5, 3) = X_35;
	SPrior(3, 6) = X_36;
	SPrior(6, 3) = X_36;
	SPrior(4, 4) = X_44;
	SPrior(4, 4) = X_44;
	SPrior(4, 5) = X_45;
	SPrior(5, 4) = X_45;
	SPrior(4, 6) = X_46;
	SPrior(6, 4) = X_46;
	SPrior(5, 5) = X_55;
	SPrior(5, 5) = X_55;
	SPrior(5, 6) = X_56;
	SPrior(6, 5) = X_56;
	SPrior(6, 6) = X_66;
	SPrior(6, 6) = X_66;
	KPrior(1) = X_146;
	KPrior(2) = X_246;
	KPrior(3) = X_346;
	KPrior(4) = X_446;
	KPrior(5) = X_456;
	KPrior(6) = X_466;
	KPrior(7) = X_4466;
	KPrior(8) = X_126;
	KPrior(9) = X_136;
	KPrior(10) = X_226;
	KPrior(11) = X_236;
	KPrior(12) = X_256;
	KPrior(13) = X_266;
	KPrior(14) = X_336;
	KPrior(15) = X_356;
	KPrior(16) = X_366;
	KPrior(17) = X_2466;
	KPrior(18) = X_3466;
	UPrior(1) = X_116;
	UPrior(2) = X_156;
	UPrior(3) = X_166;
	UPrior(4) = X_1466;
	UPrior(5) = X_2266;
	UPrior(6) = X_2366;
	UPrior(7) = X_3366;
	UPrior(8) = X_1266;
	UPrior(9) = X_1366;
	UPrior(10) = X_1166;
	MuPriorLN = 2*log(MuPrior) - 1/2*log(diag(SPrior));
	SigmaPriorLN = log(SPrior) - log(MuPrior)*E' - E*log(MuPrior)';
	InvSigmaPriorLN = inv(SigmaPriorLN);
	InvSigmaLN = w*w' / v^2 + InvSigmaPriorLN;
	SigmaPosteriorLN = inv(InvSigmaLN);
	MuPosteriorLN = SigmaPosteriorLN * (w*y / v^2 + InvSigmaPriorLN * MuPriorLN);
	m = w'*MuPriorLN;
	s = SigmaPriorLN(4, 4);
	Likelihood = -1/2*(y-m)^2/(v^2 + s) - 1/2*log(1/v^2 + 1/s) - 1/2*log(s) - log(v);
	if isnan(Likelihood) || isinf(Likelihood)
		Likelihood = -inf;
		MuPosterior = zeros(size(MuPosterior));
		SPosterior = zeros(size(SPosterior));
		SigmaPosterior = zeros(size(SigmaPosterior));
	else
	end
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
	SPosteriorZ(1, 4) = exp(MuPosteriorLN(1) + MuPosteriorLN(4) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 4)+...
		SigmaPosteriorLN(4, 1) + SigmaPosteriorLN(4, 4)));
	SPosteriorZ(4, 1) = SPosteriorZ(1, 4);
	SPosteriorZ(1, 5) = exp(MuPosteriorLN(1) + MuPosteriorLN(5) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 5)+...
		SigmaPosteriorLN(5, 1) + SigmaPosteriorLN(5, 5)));
	SPosteriorZ(5, 1) = SPosteriorZ(1, 5);
	SPosteriorZ(1, 6) = exp(MuPosteriorLN(1) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 6)+...
		SigmaPosteriorLN(6, 1) + SigmaPosteriorLN(6, 6)));
	SPosteriorZ(6, 1) = SPosteriorZ(1, 6);
	SPosteriorZ(2, 2) = exp(MuPosteriorLN(2) + MuPosteriorLN(2) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 2)+...
		SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 2)));
	SPosteriorZ(2, 2) = SPosteriorZ(2, 2);
	SPosteriorZ(2, 3) = exp(MuPosteriorLN(2) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 3)+...
		SigmaPosteriorLN(3, 2) + SigmaPosteriorLN(3, 3)));
	SPosteriorZ(3, 2) = SPosteriorZ(2, 3);
	SPosteriorZ(2, 4) = exp(MuPosteriorLN(2) + MuPosteriorLN(4) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 4)+...
		SigmaPosteriorLN(4, 2) + SigmaPosteriorLN(4, 4)));
	SPosteriorZ(4, 2) = SPosteriorZ(2, 4);
	SPosteriorZ(2, 5) = exp(MuPosteriorLN(2) + MuPosteriorLN(5) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 5)+...
		SigmaPosteriorLN(5, 2) + SigmaPosteriorLN(5, 5)));
	SPosteriorZ(5, 2) = SPosteriorZ(2, 5);
	SPosteriorZ(2, 6) = exp(MuPosteriorLN(2) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 6)+...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 6)));
	SPosteriorZ(6, 2) = SPosteriorZ(2, 6);
	SPosteriorZ(3, 3) = exp(MuPosteriorLN(3) + MuPosteriorLN(3) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3)+...
		SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3)));
	SPosteriorZ(3, 3) = SPosteriorZ(3, 3);
	SPosteriorZ(3, 4) = exp(MuPosteriorLN(3) + MuPosteriorLN(4) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 4)+...
		SigmaPosteriorLN(4, 3) + SigmaPosteriorLN(4, 4)));
	SPosteriorZ(4, 3) = SPosteriorZ(3, 4);
	SPosteriorZ(3, 5) = exp(MuPosteriorLN(3) + MuPosteriorLN(5) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 5)+...
		SigmaPosteriorLN(5, 3) + SigmaPosteriorLN(5, 5)));
	SPosteriorZ(5, 3) = SPosteriorZ(3, 5);
	SPosteriorZ(3, 6) = exp(MuPosteriorLN(3) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 6)+...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 6)));
	SPosteriorZ(6, 3) = SPosteriorZ(3, 6);
	SPosteriorZ(4, 4) = exp(MuPosteriorLN(4) + MuPosteriorLN(4) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 4)+...
		SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 4)));
	SPosteriorZ(4, 4) = SPosteriorZ(4, 4);
	SPosteriorZ(4, 5) = exp(MuPosteriorLN(4) + MuPosteriorLN(5) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 5)+...
		SigmaPosteriorLN(5, 4) + SigmaPosteriorLN(5, 5)));
	SPosteriorZ(5, 4) = SPosteriorZ(4, 5);
	SPosteriorZ(4, 6) = exp(MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6)+...
		SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));
	SPosteriorZ(6, 4) = SPosteriorZ(4, 6);
	SPosteriorZ(5, 5) = exp(MuPosteriorLN(5) + MuPosteriorLN(5) + ...
		1/2*(SigmaPosteriorLN(5, 5) + SigmaPosteriorLN(5, 5)+...
		SigmaPosteriorLN(5, 5) + SigmaPosteriorLN(5, 5)));
	SPosteriorZ(5, 5) = SPosteriorZ(5, 5);
	SPosteriorZ(5, 6) = exp(MuPosteriorLN(5) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(5, 5) + SigmaPosteriorLN(5, 6)+...
		SigmaPosteriorLN(6, 5) + SigmaPosteriorLN(6, 6)));
	SPosteriorZ(6, 5) = SPosteriorZ(5, 6);
	SPosteriorZ(6, 6) = exp(MuPosteriorLN(6) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6)+...
		SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6)));
	SPosteriorZ(6, 6) = SPosteriorZ(6, 6);

	MuPosterior=MuPosteriorZ;
	SPosterior=SPosteriorZ;

	KPosterior(1) = exp(MuPosteriorLN(1) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 4) + SigmaPosteriorLN(1, 6) +...
		SigmaPosteriorLN(4, 1) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 1) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

	KPosterior(2) = exp(MuPosteriorLN(2) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 4) + SigmaPosteriorLN(2, 6) +...
		SigmaPosteriorLN(4, 2) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

	KPosterior(3) = exp(MuPosteriorLN(3) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 4) + SigmaPosteriorLN(3, 6) +...
		SigmaPosteriorLN(4, 3) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

	KPosterior(4) = exp(MuPosteriorLN(4) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) +...
		SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

	KPosterior(5) = exp(MuPosteriorLN(4) + MuPosteriorLN(5) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 5) + SigmaPosteriorLN(4, 6) +...
		SigmaPosteriorLN(5, 4) + SigmaPosteriorLN(5, 5) + SigmaPosteriorLN(5, 6) + ...
		SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 5) + SigmaPosteriorLN(6, 6)));

	KPosterior(6) = exp(MuPosteriorLN(4) + MuPosteriorLN(6) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + SigmaPosteriorLN(4, 6) +...
		SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6) + ...
		SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6)));

	KPosterior(7) = exp(MuPosteriorLN(4) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) +...
		SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

	KPosterior(8) = exp(MuPosteriorLN(1) + MuPosteriorLN(2) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 2) + SigmaPosteriorLN(1, 6) +...
		SigmaPosteriorLN(2, 1) + SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 6) + ...
		SigmaPosteriorLN(6, 1) + SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 6)));

	KPosterior(9) = exp(MuPosteriorLN(1) + MuPosteriorLN(3) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(1, 1) + SigmaPosteriorLN(1, 3) + SigmaPosteriorLN(1, 6) +...
		SigmaPosteriorLN(3, 1) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 6) + ...
		SigmaPosteriorLN(6, 1) + SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 6)));

	KPosterior(10) = exp(MuPosteriorLN(2) + MuPosteriorLN(2) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 6) +...
		SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 6) + ...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 6)));

	KPosterior(11) = exp(MuPosteriorLN(2) + MuPosteriorLN(3) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 3) + SigmaPosteriorLN(2, 6) +...
		SigmaPosteriorLN(3, 2) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 6) + ...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 6)));

	KPosterior(12) = exp(MuPosteriorLN(2) + MuPosteriorLN(5) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 5) + SigmaPosteriorLN(2, 6) +...
		SigmaPosteriorLN(5, 2) + SigmaPosteriorLN(5, 5) + SigmaPosteriorLN(5, 6) + ...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 5) + SigmaPosteriorLN(6, 6)));

	KPosterior(13) = exp(MuPosteriorLN(2) + MuPosteriorLN(6) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 6) + SigmaPosteriorLN(2, 6) +...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6) + ...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6)));

	KPosterior(14) = exp(MuPosteriorLN(3) + MuPosteriorLN(3) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 6) +...
		SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 6) + ...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 6)));

	KPosterior(15) = exp(MuPosteriorLN(3) + MuPosteriorLN(5) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 5) + SigmaPosteriorLN(3, 6) +...
		SigmaPosteriorLN(5, 3) + SigmaPosteriorLN(5, 5) + SigmaPosteriorLN(5, 6) + ...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 5) + SigmaPosteriorLN(6, 6)));

	KPosterior(16) = exp(MuPosteriorLN(3) + MuPosteriorLN(6) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 6) + SigmaPosteriorLN(3, 6) +...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6) + ...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 6) + SigmaPosteriorLN(6, 6)));

	KPosterior(17) = exp(MuPosteriorLN(2) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(2, 2) + SigmaPosteriorLN(2, 4) + SigmaPosteriorLN(2, 6) +...
		SigmaPosteriorLN(4, 2) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 2) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

	KPosterior(18) = exp(MuPosteriorLN(3) + MuPosteriorLN(4) + MuPosteriorLN(6) + ...
		1/2*(SigmaPosteriorLN(3, 3) + SigmaPosteriorLN(3, 4) + SigmaPosteriorLN(3, 6) +...
		SigmaPosteriorLN(4, 3) + SigmaPosteriorLN(4, 4) + SigmaPosteriorLN(4, 6) + ...
		SigmaPosteriorLN(6, 3) + SigmaPosteriorLN(6, 4) + SigmaPosteriorLN(6, 6)));

