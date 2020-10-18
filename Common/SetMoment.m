function MNew = SetMoment(momentSystem, M, stateIdxVec, specIdx, value)


order = length(specIdx);
numStates = length(momentSystem);

if (isempty(stateIdxVec))
  stateIdxVec = 1;
  numStates = 0;
end


for u=1:length(stateIdxVec)
    
    stateIdx = stateIdxVec(u);
    numMoments = length(momentSystem{stateIdx}.dM);
    
    for k=1:length(momentSystem{stateIdx}.dM)
        if (momentSystem{stateIdx}.Order(k)==order)
            specIdxTmp = momentSystem{stateIdx}.SpeciesIdx{k};
            if (sum(specIdx==specIdxTmp)==order)
                idx = k;
                break;
            end
        end
    end
    
    totalIdx = numStates + (stateIdx-1)*numMoments + idx;
    M(totalIdx) = value;
    
end

MNew = M;

end