function [commonParams] = GetModelParameters(model, config)

numQParams = length(config.UnknownQIdx);
numZParams = length(config.UnknownZIdx);
numCParams = length(config.UnknownCIdx);
numHParams = length(config.UnknownHIdx);
numMParams = length(config.UnknownMIdx);

commonParams = zeros(numQParams+numZParams+numCParams+numHParams+numMParams, 1);

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


end