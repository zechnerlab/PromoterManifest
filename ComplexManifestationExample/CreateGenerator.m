function Q = CreateGenerator(config, t)

N = config.numPromoterStates;
Q = zeros(N, N);

u = GetPieceWiseConstantInput(t, config.InputParams);

for k=1:length(config.rates)
   
    from = config.rates{k}.Idx(1);
    to = config.rates{k}.Idx(2);
    
    Q(to, from) = config.rates{k}.Fun(t, u, config.rates{k}.Params);
    
end


Q = Q - diag(sum(Q));