function [discretePath, discreteTime] = SampleCTMPPath(xPath, tPath, deltaT, maxTime)


    if nargin < 4
        maxTime = inf;
    end

    % Number of species.
    M = size(xPath, 1);

    %% TODO: Sanity checks!
    minTime = min(tPath);
    maxTime = min(maxTime, max(tPath));
    
    discreteTime = minTime:deltaT:maxTime;
    N = length(discreteTime);
    
    discretePath = zeros(M, N);
    discretePath(:, 1) = xPath(:, 1);
    
    tPathIndex = 2;
    
    
    for i=2:length(discreteTime)
       
        currTime = discreteTime(i);
        
        %tPathIndex = tPathIndex + 1;
        while(currTime > tPath(tPathIndex))
           tPathIndex = tPathIndex + 1; 
        end
        
        tPathIndex = tPathIndex - 1;
        discretePath(:, i) = xPath(:, tPathIndex);
        
        
    end
end