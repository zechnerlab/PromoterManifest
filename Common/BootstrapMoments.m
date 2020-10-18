function [means, vars] = BootstrapMoments(data, N, M)

    K = length(data);
    
    for i=1:M
        rndIdx = randi(K, N, 1);
        means(i) = mean(data(rndIdx));
        vars(i) = var(data(rndIdx));
    end
    

end