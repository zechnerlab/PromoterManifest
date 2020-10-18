function QNet = CalculateFlows(P, tGrid, config)

N = config.numPromoterStates;

QNet = zeros(N, N);

for k=1:length(tGrid)-1
    
    dt = tGrid(k+1) - tGrid(k);
    Q = CreateGenerator(config, tGrid(k));
    
    PMat = repmat(P(k, :), N, 1);
    QNet = QNet + Q.*PMat*dt;
    
end