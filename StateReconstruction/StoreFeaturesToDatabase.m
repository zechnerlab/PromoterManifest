function dataB = StoreFeaturesToDatabase(dataB, features, promoter, condition)

if (isempty(fieldnames(dataB))) 
    featureNames = fieldnames(features);
    dataB.featureNames = featureNames;
    
    dataB.Data = zeros(0, length(featureNames));
    dataB.Promoter = [];
    dataB.Condition = zeros(0, 4);
end



for k=1:length(dataB.featureNames)
    featureVec(k) = real(features.(dataB.featureNames{k}));
end

dataB.Data(end+1, :) = featureVec;
dataB.Promoter{end+1} = promoter;
if (length(condition.PulseParameters)==2)
    condition.PulseParameters(3) = -1;
end
dataB.Condition(end+1, :) = [condition.Concentration, condition.PulseParameters];
end