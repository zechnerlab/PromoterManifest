function [moments, uncertainties, quantiles5, quantiles95] = BootstrapMomentsStd(dataVec, M)

L = length(dataVec);

for l=1:L
    data = dataVec{l};

    K = length(data);
    
    for i=1:M
        rndIdx = randi(K, K, 1);
        means(i) = mean(data(rndIdx));
        vars(i) = var(data(rndIdx));
    end
    
    
    moments(1, l) = mean(means);
    moments(2, l) = mean(vars);
    
    
    uncertainties(1, l) = var(means);
    uncertainties(2, l) = var(vars);
    
    covM = cov(means, vars);
    
    uncertainties(3, l) = covM(1, 2);
    
    quantiles5(1, l) = quantile(means, 0.05);
    quantiles5(2, l) = quantile(vars, 0.05);
    
    quantiles95(1, l) = quantile(means, 0.95);
    quantiles95(2, l) = quantile(vars, 0.95);
end


end