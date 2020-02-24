function logP = EvaluatePrior(model, config)

params = GetModelParameters(model, config);

numQ = length(config.UnknownQIdx);
numZ = length(config.UnknownZIdx);

q = zeros(numQ, 1);

for k=1:numQ
    idx = config.UnknownQIdx{k};
    q(k) = model.Q(idx(1), idx(2));
end

%z = zeros(numZ, 1);
z = model.Z(config.UnknownZIdx);
z = z(:);
alpha = model.alpha(:);
beta = model.beta(:);

if (strcmp(config.QPrior, 'Laplace'))
    
    if (isempty(z))
        z = 0;
        beta = 1;
    end
    
    logP = -2/alpha*sum(abs(q)) - log(alpha) - 2/beta*sum(abs(z)) - log(beta);
    %logP = -2/alpha*sum(abs(1./diag(model.Q))) - log(alpha) - 2/beta*sum(abs(model.Z)) - log(beta);
elseif (strcmp(config.QPrior, 'LaplaceNull'))
    
    n = null(model.Q);
    n = n/sum(n);
    
    logP = -2/alpha*sum(n) - log(alpha) - 2/beta*sum(abs(model.Z)) - log(beta);
elseif (strcmp(config.QPrior, 'GammaAll'))
    
    allParams = [q; z];
    
    logP = sum(log(gampdf(allParams, alpha, 1/beta)));
elseif (strcmp(config.QPrior, 'LaplaceGroup'))
    
    colSum = sum(abs(model.Q - diag(diag(model.Q))), 2);
    
    logP = sum(-2./alpha.*colSum - log(alpha)) + sum(-2./beta.*abs(z) - log(beta));
elseif (strcmp(config.QPrior, 'LaplaceAll'))
    
    %alpha = params(numQ+numZ+1:numQ+numZ+numQ);
    %beta = params(end-numZ+1:end);
    
    if (isempty(z))
       z = 0;
       beta = 1;
    end
    
    if (sum([q;z]<=0)>0)
        logP = -inf;
    else
        logP = sum(-2./alpha.*abs(q) - log(alpha)) + sum(-2./beta.*abs(z) - log(beta));
    end
elseif (strcmp(config.QPrior, 'GaussAll'))
    
    alpha = params(numQ+numZ+1:numQ+numZ+numQ);
    beta = params(end-numZ+1:end);
    
    logP = sum(-1/2*power(q./alpha, 2) - log(alpha)) + sum(-1/2*power(z./beta, 2) - log(beta));% - sum(log(beta)) - sum(log(alpha));
 elseif (strcmp(config.QPrior, 'None'))
    
   logP = 0;
   
 elseif (strcmp(config.QPrior, 'Exponential'))
    
   logP = sum(-1./alpha.*abs(q) - log(alpha)) + sum(-1./beta.*abs(z) - log(beta));

      
 elseif (strcmp(config.QPrior, 'Gamma'))
    
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


if (~isempty(config.UnknownMIdx))
   a = 8;
   meanV = 2000;
   b = a/meanV;
   logP = logP + log(gampdf(model.MeasurementSigma, a, 1/b));
end


