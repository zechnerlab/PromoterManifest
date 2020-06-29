function kineticPart = ModifyZeroOrderFluxes(kineticPart, X)

if (~isequaln(kineticPart, sym(0)))
    
    zeroOrder = 0;
    zoCount = 0;
    for u=1:length(X)
        
        [polCoeff, T] = coeffs(kineticPart, X(u));
        
        eVec = [];
        for n=1:length(T)
            eVec(n) = isequal(T(n), sym(1));
        end
        
        if sum(eVec)>0
            if ~isequal(polCoeff(end), kineticPart)
                zeroOrder = 1;
                zeroOrderCoeff = polCoeff(end);
                break;
            else
                zoCount = zoCount + 1;
            end
        end
        
    end
    
    if zoCount==length(X)
       zeroOrder=1; 
       zeroOrderCoeff = polCoeff(end);
    end
    
    if (zeroOrder==1)
        kineticPart = simplify(kineticPart - zeroOrderCoeff + zeroOrderCoeff*sym('PM'));
    end
end

