function logP = EvaluatePrior(model, config)


numQ = length(config.UnknownQIdx);
q = zeros(numQ, 1);

for k=1:numQ
    idx = config.UnknownQIdx{k};
    q(k) = model.Q(idx(1), idx(2));
end

z = model.Z(config.UnknownZIdx);
z = z(:);


%% Evalulate priors for switching rates q_ij and transcription rates z_1, z_2
      
if (strcmp(config.QPrior, 'Gamma'))
   %hard-coded prior with alpha=1 and beta=30, corresponding to an average
   %switching and transcription rate of 1/30 s^-1. With alpha=1, the prior
   %is not strongly informative.
   alpha = 1;
   beta = 30;
   logP = sum(log(gampdf(q, alpha, 1/beta))) + sum(log(gampdf(z, alpha, 1/beta)));
   
end


%% evaluate priors for rate constants

logCPrior = [];
for k=1:length(config.UnknownCIdx)
   logCPrior(k) = log(gampdf(model.c(config.UnknownCIdx(k)), config.CPriors.A(k), 1/config.CPriors.B(k))); 
end

logP = logP + sum(logCPrior);



