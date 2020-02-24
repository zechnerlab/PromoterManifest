function [h, g] = CalculatePropensities(x, c, Pre)

    % M = Number of Species, N = Number Jumps.
    [M, N] = size(x);

    % Number of Reactions.
    K = length(c);
    
    % Initialize propensity vector.
    h = zeros(K, N);
    g = zeros(K, N);
    
   
    for k=1:K
        
       % For some known and simple reactions, we can compute the propensities
       % efficiently without calculating binomial coefficients. 
        
       %pre = Pre(k, :);
       
       %%%%%!!!! Commented out sorting, so complicated reactions are no
       %%%%%more supported...was necessary to achieve some significant
       %%%%%speedup!
       %% Sort multiplicities of the reactants.
       %[sortedPre, sortedPreIndex] = sort(pre, 'descend');
       
       % Find out involved reactants (multiplicity > 0).
       reactantIndex = logical(Pre(k, :));% > 0;
       multiplicities = Pre(k, reactantIndex);
       numReactants = length(multiplicities);
       
       
       if (numReactants == 0)
           g(k, :) = 1;
           continue;
       elseif (numReactants == 1)
           
           numMolecules = x(reactantIndex, :);
           
           if (multiplicities == 1)
               g(k, :) = numMolecules;
               continue;
           elseif (multiplicities == 2)
               g(k, :) = 0.5 * numMolecules .* (numMolecules - 1);
               continue;
           elseif (multiplicities == 3)
               g(k, :) = 1/6 * numMolecules .* (numMolecules - 1)...
                   .* (numMolecules - 2);
               continue;
           end
       elseif (numReactants == 2)
           [numMolecules] = x(reactantIndex, :);
           
           if (multiplicities(1) == 1 && multiplicities(2) == 1)
               g(k, :) = numMolecules(1, :) .* numMolecules(2, :);
               continue;
           elseif (multiplicities(1) == 2 && multiplicities(2) == 1)
               g(k, :) = 0.5 * numMolecules(1, :) .* (numMolecules(1, :) - 1)...
                   .* numMoleculesB;
               continue;
           end
       end
       
       
       %% For all other reactions...
       for i=1:N
            
            stochProd = 1;
           
            % For each species...
            for l=1:M

                % check if there are at least as many molecules as required
                % for the reaction to take place.
                if (x(l, i) >= Pre(k, l))

                    % ...calculate factor depending on the multiplicity of
                    % the current species. This is give by the left-hand 
                    % side of the reaction. The rounding is necessary if 
                    % the states have been reconstructed using real 
                    % numbers!
                    s = x(l, i);
                    stochProd = stochProd * nchoosek(round(x(l, i)), ...
                        Pre(k, l));
                elseif (x(l, i) == 0 && Pre(k, l) > 1)
                    % ...set the propensity to zero and abort if muliplicity is
                    % too small for the reaction to take place.
                    stochProd = 0;
                    break;
                end
            end

            % Multiply by stochastic rate constant.
            g(k, i) = stochProd;
        end
       
    end
  
    h = g .* c(:, ones(1, N));
end