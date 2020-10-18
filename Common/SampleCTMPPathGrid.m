function [discretePath, discreteTime] = SampleCTMPPathGrid(xPath, tPath, grid)


    if nargin < 4
        maxTime = inf;
    end

    % Number of species.
    M = size(xPath, 1);

    %% TODO: Sanity checks!
    
    discreteTime = grid;
    N = length(discreteTime);
    
    discretePath = zeros(M, N);
    
    tPathIndex = 1;
    
    
    for i=1:length(discreteTime)
       
        currTime = discreteTime(i);
        
        %tPathIndex = tPathIndex + 1;
        
        while(currTime > tPath(tPathIndex))
           tPathIndex = tPathIndex + 1; 
           
           if (tPathIndex > length(tPath))
               break;
           end
        end
        
        discretePath(:, i) = xPath(:, tPathIndex);
        
        
    end
end