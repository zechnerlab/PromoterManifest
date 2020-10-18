%% SimulateGillespieTimeVaryingSync - Simulate sample path of a reaction 
% system with a time-varying propensity function. The time-varying
% propensity function is updated (i.e., synchronized at given discrete time
% points syncTime. X0 is a row-vector containing the initial state of the
% reaction system and c is a K-by-1 matrix containing the reaction rate constants
% corresponding to the K reaction channels. Pre and Post are K-by-L
% matrices such that the stoichimetry matrix is given by S=Post-Pre.
% and Post are given by the stoichiometry. N defines the maximum number of
% reactions and maxTime is the final simulation time. The reaction with
% index tvRateIndex is associated with a time-varying rate constant which
% is specified by a function handle funHandle and corresponding parameters
% params. Additionally, the start time of the simulation can be specified
% by startTime.
%
% [X, t, r, H, G, cVarying] = SimulateGillespieTimeVaryingSync(X0,
%       c, Pre, Post, N, maxTime, 
%       tvRateIndex, funHandle, 
%       params, syncTime, startTime)

function [X, t, r, H, G, cVarying] = SimulateGillespieTimeVaryingSync(X0, c, Pre,...
                    Post, N, maxTime, tvRateIndex, funHandle, params, ...
                    syncTime, startTime)

             
    if (nargin < 11)
       startTime = 0; 
    end


    % Number of species.
    M = length(X0);
    
    % Number of reactions.
    K = length(c);

    X = zeros(M, N+1);
    X(:, 1) = X0;
    
    % Calculate stoichiometric matrix.
    S = Post - Pre;
    
    t = zeros(1, N+1);
    t(1) = startTime;
    
    r = zeros(1, N+1);
    
    
    H = zeros(K, N);
    G = zeros(K, N);
    
    cVarying = zeros(length(tvRateIndex), N);
    
    % if T0 is not zero, initialize syncCounter
    for i=1:length(syncTime)
       if (syncTime(i) >= startTime)
          syncCounter = i;
          break;
       end
    end
    
    c(tvRateIndex) = funHandle(startTime, params);
    
    for i=1:N
        
        if (t(i) >= maxTime)
            % Cut away last reaction
            X = [X(:, 1:i-1), X(:, i-1)];
            t = [t(1:i-1), maxTime];
            H = H(:, 1:i-1);
            G = G(:, 1:i-1);
            r = r(1:i-1);
            cVarying = cVarying(:, 1:i-1);
            break;
        end
        

       % Calculate propensities for the current state vector.
       [h, g] = CalculatePropensitiesFast(X(:, i), c, Pre);  
        
       cVarying(:, i) = c(tvRateIndex);
       
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
           cVarying = cVarying(:, 1:i);
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
           
           if (t(i+1) > syncTime(syncCounter))
              c(tvRateIndex) = funHandle(syncTime(syncCounter), params);
              t(i+1) = syncTime(syncCounter);
              syncCounter = syncCounter + 1;
              X(:, i+1) = X(:, i);
           end
           
       end
       
       X(X(:, i+1)<0, i+1) = 0;
       
       r(i) = reactionIndex;
    end

end