function model = SetModelParameters(model, config, commonParams)

    numQParams = length(config.UnknownQIdx);
    numZParams = length(config.UnknownZIdx);
    numCParams = length(config.UnknownCIdx);
    numHParams = length(config.UnknownHIdx);
    numMParams = length(config.UnknownMIdx);
    numIParams = 0;%1;%size(model.Q, 1)-2;
    numXParams = 0;%2;
    
    if (config.InferHyper==1)
        if (strcmp(config.QPrior, 'LaplaceAll'))
            numAParams = numQParams;
            numBParams = numZParams;
        elseif (strcmp(config.QPrior, 'LaplaceGroup'))
            numAParams = length(model.P0);
            numBParams = numZParams;
        else
            numAParams = min(numQParams, 1);
            numBParams = min(numZParams, 1);
        end
    else
        numAParams = 0;
        numBParams = 0;
    end

    QParams = commonParams(1:numQParams);
    ZParams = commonParams(numQParams+1:numQParams+numZParams);
    CParams = commonParams(numQParams+numZParams+1:numQParams+numZParams+numCParams);
    HParams = commonParams(numQParams+numZParams+numCParams+1:numQParams+numZParams+numCParams+numHParams);
    MParams = commonParams(numQParams+numZParams+numCParams+numHParams+1:numQParams+numZParams+numCParams+numHParams+numMParams);
    AParams = commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+1:numQParams+numZParams+numCParams+numHParams+numMParams+numAParams);
    BParams = commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+1:numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams);
    IParams = commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+1:numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+numIParams);
    XParams = commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+numIParams+1:numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+numIParams+numXParams);
    
    for l=1:numQParams
       idx = config.UnknownQIdx{l};
       model.Q(idx(1), idx(2)) = QParams(l);
    end
    
    model.Q = model.Q - diag(sum(model.Q, 1));
    
    model.Z(config.UnknownZIdx) = ZParams;
    model.c(config.UnknownCIdx) = CParams;
    model.H(config.UnknownHIdx) = HParams;
    model.MeasurementSigma(config.UnknownMIdx) = MParams;
    
    if (config.InferHyper==1)
        model.alpha = AParams;
        model.beta = BParams;
    end
    
    %N = null(model.Q);
    %model.P0 = N/sum(N);
    %model.P0(2) = IParams;
    %model.P0(1) = 0;
    %model.P0(1) = 1-sum(model.P0);
    %model.P0(end) = 0;
    
    muExt = model.H(1);
    sigmaExt_squared = model.H(2);
    X0 = model.M0;
    numStates = size(model.Q, 1);
    odeInfos = model.odeInfos;
    
    X0(1:numStates) = model.P0;%1/numStates;
    
%     X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], 3, muExt*1/numStates); %set E[Z]
%     X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [3 3], (muExt^2+sigmaExt_squared)*1/numStates); %set E[ZZ]
%     X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [2], model.muMeas0*1/numStates);
%     X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [2 2], (model.muMeas0^2 + model.varMeas0)*1/numStates);
%     X0 = SetMoment(odeInfos.infos.MomentSystem, X0, [1:numStates], [2 3], model.muMeas0*muExt*1/numStates);
    
    %model.InitialMeans = XParams(1);
    %model.InitialVars = XParams(2);

    means = [model.InitialMeans; model.muMeas0; muExt];
    vars = [model.InitialVars; model.varMeas0; sigmaExt_squared];
    
    X0 = SetInitialMoments(odeInfos.infos, X0, numStates, means, vars);
    
%     for k=1:numStates
%         X0 = SetMoment(odeInfos.infos.MomentSystem, X0, k, 4, muExt*X0(k)); %set E[Z]
%         X0 = SetMoment(odeInfos.infos.MomentSystem, X0, k, [4 4], (muExt^2+sigmaExt_squared)*X0(k)); %set E[ZZ]
%         X0 = SetMoment(odeInfos.infos.MomentSystem, X0, k, [3], model.muMeas0*X0(k));
%         X0 = SetMoment(odeInfos.infos.MomentSystem, X0, k, [3 3], (model.muMeas0^2 + model.varMeas0)*X0(k));
%         X0 = SetMoment(odeInfos.infos.MomentSystem, X0, k, [3 4], model.muMeas0*muExt*X0(k));
%     end
    
    model.M0 = X0;

    
    %if (config.InferHyper == 1)
    %   model.alpha = params(numQParams+numZParams+1:numQParams+numZParams+numQParams);
    %   model.beta = params(end-numZParams+1:end);
    %end

    
end