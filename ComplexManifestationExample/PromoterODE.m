function dX = PromoterODE(t, X, config)

N = config.numPromoterStates;
Z = config.Z;


P = X(1:N);
Moments = X(N+1:end);
m = Moments(1:N);
s = Moments(N+1:end);

Q = CreateGenerator(config, t);

dP = Q*P;


dm = diag(Z)*P + Q*m;
ds = diag(Z)*P + 2*diag(Z)*m + Q*s;

dX = [dP(:); [dm(:); ds(:)]];

end