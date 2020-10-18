function model = SetModelKineticParamsOnly(model, config, commonParams)

    numQParams = length(config.UnknownQIdx);
    numZParams = length(config.UnknownZIdx);
    numCParams = length(config.UnknownCIdx);
    numHParams = length(config.UnknownHIdx);
    numMParams = length(config.UnknownMIdx);


    QParams = commonParams(1:numQParams);
    ZParams = commonParams(numQParams+1:numQParams+numZParams);
    CParams = commonParams(numQParams+numZParams+1:numQParams+numZParams+numCParams);
    HParams = commonParams(numQParams+numZParams+numCParams+1:numQParams+numZParams+numCParams+numHParams);
    MParams = commonParams(numQParams+numZParams+numCParams+numHParams+1:numQParams+numZParams+numCParams+numHParams+numMParams);
   
    
    for l=1:numQParams
       idx = config.UnknownQIdx{l};
       model.Q(idx(1), idx(2)) = QParams(l);
    end
    
    model.Q = model.Q - diag(sum(model.Q, 1));
    
    model.Z(config.UnknownZIdx) = ZParams;
    %model.c(config.UnknownCIdx) = CParams;
    %model.H(config.UnknownHIdx) = HParams;
    %model.MeasurementSigma(config.UnknownMIdx) = MParams;
    

    
end