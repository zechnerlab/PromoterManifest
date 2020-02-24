function [commonParams] = GetModelParameters(model, config)

numQParams = length(config.UnknownQIdx);
numZParams = length(config.UnknownZIdx);
numCParams = length(config.UnknownCIdx);
numHParams = length(config.UnknownHIdx);
numMParams = length(config.UnknownMIdx);
numIParams = 0;%size(model.Q, 1)-2;
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

commonParams = zeros(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+numIParams+numXParams, 1);

for l=1:numQParams
    idx = config.UnknownQIdx{l};
    commonParams(l) = model.Q(idx(1), idx(2));
end

for l=1:numZParams
    commonParams(numQParams+l) = model.Z(config.UnknownZIdx(l));
end

for l=1:numCParams
    commonParams(numQParams+numZParams+l) = model.c(config.UnknownCIdx(l), 1);
end

for l=1:numHParams
    commonParams(numQParams+numZParams+numCParams+l) = model.H(config.UnknownHIdx(l));
end

for l=1:numMParams
    commonParams(numQParams+numZParams+numCParams+numHParams+l) = model.MeasurementSigma(config.UnknownMIdx(l));
end


for l=1:numAParams
    commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+l) = model.alpha(l);
end

for l=1:numBParams
     commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+l) = model.beta(l);
end

for l=1:numIParams
     commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+l) = model.P0(length(model.P0)-l);
end


%commonParams(numQParams+numZParams+numCParams+numHParams+numMParams+numAParams+numBParams+numIParams+1:end) = [model.InitialMeans; model.InitialVars];



%if (config.InferHyper==1)
%   commonParams = [commonParams; model.alpha];
%   commonParams = [commonParams; model.beta];
%end

end