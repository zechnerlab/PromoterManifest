function [Pre, Post, c, X0] = CreatePromoterModel(Q)


    numStates = length(Q);        
  	X0 = zeros(numStates, 1);
    X0(1) = 1;
    
    count = 1;
    for k=1:numStates
        for i=1:numStates
            if (k==i)
                continue;
            end
            Pre(count, k) = 1;
            Post(count, i) = 1;
            
            c(count) = Q(i, k);
            count = count + 1;
        end
    end
         
    c = c(:);

end