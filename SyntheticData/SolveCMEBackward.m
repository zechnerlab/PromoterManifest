function [mRNAVec, PmRNA, PAll] = SolveCMEBackward(P0, c, M, timePoints)


% enumerate states:


gene(1:M) = 0;
gene(M+1:2*M) = 1;
mRNA(1:M) = 0:M-1;
mRNA(M+1:2*M) = 0:M-1;

mRNAVec = mRNA(1:M)';
mRNAVec = [mRNAVec; mRNAVec];

counter = 1;

%times = union(timePoints, inputParams.inputTimes);
times = timePoints;

prevT = 0;

cTmp = c;
for l=1:length(times)
    
    t = times(l);
    
    %cTmp(inputParams.inputRateIndex) = GetPieceWiseConstantInput(t, inputParams);
    geneOn = cTmp(1);
    geneOff = cTmp(2);
    mRNASynthesis = cTmp(3);
    mRNADegradation = cTmp(4);
    Q = zeros(2*M, 2*M);
    
    deltaT = t - prevT;

    
    for k=1:2*M
        g = gene(k);
        m = mRNA(k);
        
        
        if (m==0)
            if (g==1)
                Q(k, k+1) = mRNADegradation * (m+1);
                Q(k, k-M) = geneOn;
                Q(k, k) = -mRNADegradation * m - geneOff - mRNASynthesis;
                
            else
                Q(k, k+1) = mRNADegradation * (m+1);
                Q(k, k+M) = geneOff;
                Q(k, k) = -mRNADegradation * m - geneOn;
            end
        elseif (m==M-1)
            if (g==1)
                Q(k, k-1) = mRNASynthesis;
                Q(k, k-M) = geneOn;
                Q(k, k) = -mRNADegradation * m - geneOff;
                
            else
                Q(k, k+M) = geneOff;
                Q(k, k) = -mRNADegradation* m - geneOff;
            end
        else
            if (g==1)
                Q(k, k-1) = mRNASynthesis;
                Q(k, k+1) = mRNADegradation * (m+1);
                Q(k, k-M) = geneOn;
                Q(k, k) = - mRNADegradation * m - geneOff - mRNASynthesis;
                
            else
                Q(k, k+1) = mRNADegradation * (m+1);
                Q(k, k+M) = geneOff;
                Q(k, k) = - mRNADegradation * m - geneOn;
            end
        end
        
    end
%     
    Q = Q';
    Q = Q - diag(diag(Q));
    Q = Q - diag(sum(Q, 1));
   
    P = expm(Q*deltaT) * P0;
    
    P0 = P;
    
    if (sum(times(l) == timePoints)==1)
        PmRNA{counter} = P(1:M) + P(M+1:end);
        PmRNA{counter} = PmRNA{counter} / sum(PmRNA{counter});
        PAll{counter} = P;
        counter = counter + 1;
    end
    prevT = t;
end




end