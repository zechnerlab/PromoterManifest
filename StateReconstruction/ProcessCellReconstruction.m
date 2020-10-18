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
    
    
    %For each particle, store the corresponding statistics.
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
        
        %Promoter states 2 and 3 are considered active only if they have a
        %transcription rate which is at least x% of the maximum
        %transcription rate of that promoter. The variable maxZ is
        %extracted from the preceeding model calibration and handed over to
        %this function. A threshold value of 10-20% worked well.
        activeState = find(particles{k}.Model.Z>0.2*maxZ);
        
        
        %count number of switches from each promoter state to each other
        %promoter state. This is not used for the subsequent analyses but
        %may still carry interesting information.
        for j=1:6
            numRSwitches(k, j) = sum(particles{k}.r==j);
        end
        
        %Calculate the time spent in each of the three states.
        dT = diff(particles{k}.t);
        for j=1:3
           tauState(k, j) = sum(dT(find(particles{k}.stateIdx(1:end-1)==j))); 
        end
        
        
        %Calculate transcriptional output, which is just the time-integral of
        %the transcription rate.
        trOutput(k) = tauState(k, :)*particles{k}.Model.Z';
        
        
        %Calculate when the promoter was *in each* of the active states.
        S = [];
        for l=1:length(activeState)
           S(l, :) =  particles{k}.stateIdx == activeState(l);
        end
        
        %Determine when the promoter was in *any* active
        %state.
        responders = sum(S, 1)>0;
        
        %Count the number of switches between an inactive and an active
        %promoter state.
        numSwitches(k) = sum(and(responders(1:end-1)==0, responders(2:end)==1));
        
        %Find the time index of when the promoter was in an active state.
        activeIdx = find(responders);
        
        
        %The first of these indices corresponds to the time when the
        %promoter activated for the first time (i.e., 'timeToActivate').
        if (isempty(activeIdx))
            timeToActive(k) = nan;
        else
            timeToActive(k) = particles{k}.t(activeIdx(1));
        end
        
        %Calculate the total time the promoter was in *any* active state.
        deltaT = [diff(particles{k}.t), 0];
        totalTimeActive(k) = sum(deltaT(activeIdx));
    end
    
    
    %Here, the statistics are average over the different particles. The
    %particle weights are first normalized and then a wighted average is
    %calculate. Note: the way the SMC method operates now, the particles
    %are redrawn at the very end of the iteration such that this step is
    %not necessary. However, it may be necessary when the algorithm is
    %operated differently in possible future versions of the code.
    W = W/sum(W);
    
    
    %In order to suppress noise in the response detection, only particles which
    %resided in an active state for at least 2 minutes are considered
    %'responders'. Technically, the average below is a conditional expectation, where 
    %we condition on the cell being a responder.
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
    
    %Determine most likely particle, which serves as an approximate maximum
    %a posterior reconstruction of the promoter dynamics.
    [~, MAPIdx] = max(W);
    
    %Calculate how many distinct particules from the first iteration
    %survived until the last iteration. This allows to monitor whether the
    %SMC method degenerates.
    numUniqueParticles = length(unique(Idx));
    