function [Mu, Sigma, pActivated, p2, p3, meanTauS, meanTauA, meanNumSwitches, meanTr, varTr, meanTrOutput, meanRSwitches, meanTauState, numUniqueParticles, MAPIdx] = ProcessCellReconstruction(model, particles, tGrid, maxZ)

    ZMat = [];
    timeToActive = [];
    totalTimeActive = [];
    p2 = [];
    p3 = [];
    rnaMuMat = [];
    proteinMuMat = [];
    zMuMat = [];
    rnaSMat = [];
    proteinSMat = [];
    zSMat = [];
    meanTauS = [];
    meanTauA = [];
    meanRSwitches = [];
    meanTr = [];
    meanNumSwitches = [];
    meanTrOutput = [];
    
    M = length(particles);
    W = zeros(M, 1);
    
    for k=1:M
        
        Idx(k) = particles{k}.Idx;
        W(k) = particles{k}.W;
        ZMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.Z, particles{k}.t, tGrid);
        P2(k, :) = SampleCTMPPathGrid_mex(double(particles{k}.stateIdx==2), particles{k}.t, tGrid);
        P3(k, :) = SampleCTMPPathGrid_mex(double(particles{k}.stateIdx==3), particles{k}.t, tGrid);
        rnaMuMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.PosteriorStats(1, :), particles{k}.PosteriorTime, tGrid);
        proteinMuMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.PosteriorStats(2, :), particles{k}.PosteriorTime, tGrid);
        zMuMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.PosteriorStats(3, :), particles{k}.PosteriorTime, tGrid);
        
        rnaSMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.PosteriorStats(4, :), particles{k}.PosteriorTime, tGrid);
        proteinSMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.PosteriorStats(7, :), particles{k}.PosteriorTime, tGrid);
        zSMat(k, :) = SampleCTMPPathGrid_mex(particles{k}.PosteriorStats(9, :), particles{k}.PosteriorTime, tGrid);
        
        %activeState = find(and(particles{k}.Model.Z>0.2*max(particles{k}.Model.Z), particles{k}.Model.Z>0.02));
        
        activeState = find(particles{k}.Model.Z>0.2*maxZ);
        
        for j=1:6
            numRSwitches(k, j) = sum(particles{k}.r==j);
        end
        
        dT = diff(particles{k}.t);
        for j=1:3
           tauState(k, j) = sum(dT(find(particles{k}.stateIdx(1:end-1)==j))); 
        end
        
        trOutput(k) = tauState(k, :)*particles{k}.Model.Z';
        
        S = [];
        
        for l=1:length(activeState)
           S(l, :) =  particles{k}.stateIdx == activeState(l);
        end
        
        responders = sum(S, 1)>0;
        numSwitches(k) = sum(and(responders(1:end-1)==0, responders(2:end)==1));
        activeIdx = find(responders);
        if (isempty(activeIdx))
            timeToActive(k) = nan;
        else
            timeToActive(k) = particles{k}.t(activeIdx(1));
        end
        
        deltaT = [diff(particles{k}.t), 0];
        
        totalTimeActive(k) = sum(deltaT(activeIdx));
    end
    
    W = W/sum(W);
    
    validIdx = ~isnan(timeToActive);
    validIdx = totalTimeActive > 2*60;
    meanTauS = W(validIdx)'*timeToActive(validIdx)' / sum(W(validIdx));
    meanTauA = W'*totalTimeActive';
    pActivated = W'*validIdx';
    
    meanTrOutput = mean(trOutput);
    
    p2 = W'*P2;
    p3 = W'*P3;
    
    meanTr = W'*ZMat;
    varTr = W'*ZMat.^2 - meanTr.^2;
    
    Mu(1, :)= W'*rnaMuMat;
    Mu(2, :) = W'*proteinMuMat;
    Mu(3, :) = W'*zMuMat;
    
    Sigma(1, :) = W'*rnaSMat - Mu(1, :).^2;
    Sigma(2, :) = W'*proteinSMat - Mu(2, :).^2;
    Sigma(3, :) = W'*zSMat - Mu(3, :).^2;
    
    meanNumSwitches = mean(numSwitches);
    meanTr = mean(ZMat);
    
    meanRSwitches = mean(numRSwitches);
    meanTauState = mean(tauState);
    
    [~, MAPIdx] = max(W);
    numUniqueParticles = length(unique(Idx));
    
    
    %plot(tGrid, meanTr); hold on;
    %drawnow