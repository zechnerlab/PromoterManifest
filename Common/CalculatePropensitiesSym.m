function [g, h] = CalculatePropensitiesSym(x, Pre, c)

    [numReactions, numSpecies] = size(Pre);
    
    for k=1:numReactions
        for i=1:numSpecies
            product(i) = sym(x{i}).^Pre(k, i);
        end
        
        g(k) = prod(product);
        h(k) = sym(c{k})*g(k);
    end


end