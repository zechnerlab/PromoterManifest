function [X, t, r, H, G] = SimulateGillespieRL(X0, c, Pre, Post, N, ...
    maxTime, RLIdx, RLFun, RLFunParams)



    % Number of species.
    M = length(X0);
    
    % Number of reactions.
    K = length(c);

    X = zeros(M, N+1);
    X(:, 1) = X0;
    
    % Calculate stoichiometric matrix.
    S = Post - Pre;
    
    t = zeros(1, N+1);
    t(1) = 0;
    
    r = zeros(1, N+1);
    
    
    H = zeros(K, N);
    G = zeros(K, N);
    
    for i=1:N
        
        if (t(i) >= maxTime)
            % Cut away last reaction
            X = [X(:, 1:i-1), X(:, i-1)];
            t = [t(1:i-1), maxTime];
            H = H(:, 1:i-1);
            G = G(:, 1:i-1);
            r = r(1:i-1);
            break;
        end
        
        
       % Calculate propensities for the current state vector.
       [h, g] = CalculatePropensitiesFast(X(:, i), c, Pre);  
       
       for k=1:length(RLIdx)
           fun = RLFun{k};
           params = RLFunParams{k};
           h(RLIdx(k)) = fun(t(i), X(:, i), params);
       end
       
       
       H(:, i) = h;
       G(:, i) = g;
       
       % Calculate total propensity h0. Time to the next reaction is
       % distributed as deltaT ~ exp(t/h0).
       h0 = sum(h);
       
       
       % No reaction possible!!
       if (h0 <= 0 || isnan(h0) || isinf(h0))
           
           X = [X(:, 1:i), X(:, i)];
           t = [t(1:i), maxTime];
           H = H(:, 1:i);
           G = G(:, 1:i);
           r = r(1:i);
           
           break;
       else           
           % Calculate probability for each reaction.
           probs = h./h0;
           
           % Draw next reaction.
           reactionIndex = DrawSample(1, probs);

           % Update state vector according to the stoichiometric matrix.
           X(:, i+1) = X(:, i) + S(reactionIndex, :)';
          
           % Increase time by deltaT.
           deltaT = exprnd(1/h0);
           t(i+1) = t(i) + deltaT;
           
       end
       
       X(X(:, i+1)<0, i+1) = 0;
       
       r(i) = reactionIndex;
    end

end