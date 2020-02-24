function [odeFun, infos] = GenerateMomentsInputLN(modelIn)


X0 = modelIn.X0;
c = modelIn.c;
Pre = modelIn.Pre;
Post = modelIn.Post;
syms u;

numSpecies = length(X0);
model.Pre = Pre;
model.Post = Post;
model.X = sym('X_', [length(X0), 1]);
model.X0 = X0;
model.InputRateIdx = 1;
if (~isfield(modelIn, 'ObservedSpeciesIdx'))
    model.ObservedSpeciesIdx = 2;
else
    model.ObservedSpeciesIdx = modelIn.ObservedSpeciesIdx;
end

model.c = sym('c', [length(c), 1]);
syms Pm;
for i=1:length(model.c)
    %if (sum(Pre(i, :)>0))
    for k=1:length(model.X0)
        G(i, k) = nchoosek(model.X(k), Pre(i, k));
    end
    %else
    %G(i, 1) = Pm; %important for conditional moments (multiply with marginal)
    %only required for zero order reactions.
    %end
end

model.G = prod(G, 2);
model.H = model.c.*model.G;

mCount = 0;
maxOrder = 4;

allMoments = {};

for l=1:maxOrder
    
    PermMat = npermutek(1:numSpecies, l);
    PermMat = sort(PermMat, 2);
    PermMat = unique(PermMat, 'rows');
    for u=1:size(PermMat, 1)
        
        mStr = '';
        for i=1:l
            mStr = sprintf('%s%d', mStr, PermMat(u, i));
            mIdx(i) = PermMat(u, i);
        end
        
        mCount = mCount + 1;
        g = prod(model.X(PermMat(u, :)));%model.X(1)^2;
        S = Post - Pre;
        XCondStr = sprintf('X_');
        XCond = sym(XCondStr, [length(model.X0), 1]);%substr(model.X, 'X', XCondStr);
        
        switchRateStr = sprintf('c%d', model.InputRateIdx);
        H = subs(model.H, switchRateStr, 'Z');
        
        % compute part that comes from the kinetic model (standard momments
        % more or less)
        for k=1:length(c)
            rOrder = sum(Pre(k, :));
            gPlus = subs(g, model.X, model.X+S(k, :)');
            kineticPart = simplify(gPlus*H(k) - g*H(k));
            kineticPart = subs(kineticPart, model.X, XCond);
            kineticParts(k) = kineticPart;
        end
        kineticPartTot = expand(sum(kineticParts));
        
        % compute part that is due to switching (cross terms with other
        % conditional moments)
        
        totalMomentDynamics = kineticPartTot;
        moments{1}.dM(mCount) = totalMomentDynamics;
        moments{1}.Names{mCount} = sprintf('X_%s', mStr);
        moments{1}.SpeciesIdx{mCount} = mIdx;
        moments{1}.Order(mCount) = length(mIdx);
        allMoments{end+1} = sprintf('X_%s', mStr);
        
        
        %moments{1}.dM = SubstituteMoment(XCond(:, i), moments{1}.dM, maxOrder);
        %     totalMomentDynamics = SubstituteMoment(XCond(:, i), totalMomentDynamics, maxOrder);
    end
end

warning off;
moments{1}.dM = SubstituteMoment(XCond, moments{1}.dM, maxOrder);



%% plug in independence assumptions + closures
% for k=1:numStates
%     for i=1:numStates
%         moments{k}.dM = ApplyIndependence(moments{k}.dM, 1, 3);
%     end
% end

%% remove redundant equations

% base moments = moments of order 1 and 2


allMomentIdx = EnumerateMomentIdx(numSpecies, 1:2);

u = 1;
momentIdx = EnumerateMomentIdx(numSpecies, 1:2);
numOrderTwoMoments = length(momentIdx);
numOrderThreeMoments = length(moments{u}.Names);

dMTmp = moments{u}.dM(1:numOrderTwoMoments);
NamesTmp= moments{u}.Names(1:numOrderTwoMoments);
IdxTmp = moments{u}.SpeciesIdx(1:numOrderTwoMoments);
OrderTmp = moments{u}.Order(1:numOrderTwoMoments);

NamesToTry = moments{u}.Names(numOrderTwoMoments:end);

l = numOrderTwoMoments;

% Has to be repeated at least K times wheres K is the number of third order
% moments. Not the most effective method but works.

for j=1:(numOrderThreeMoments-numOrderTwoMoments)
    for l=(numOrderTwoMoments+1):numOrderThreeMoments
        candidate = moments{u}.Names{l};
        if (~isempty(strfind(char(dMTmp), candidate)) && isempty(strfind(cell2mat(NamesTmp), candidate)))
            dMTmp(end+1) = moments{u}.dM(l);
            NamesTmp{end+1} = moments{u}.Names{l};
            IdxTmp{end+1} = moments{u}.SpeciesIdx{l};
            OrderTmp(end+1) = moments{u}.Order(l);
        end
    end
end

moments{u}.dM = dMTmp;
moments{u}.Names = NamesTmp;
moments{u}.SpeciesIdx = IdxTmp;
moments{u}.Order = OrderTmp;



%% Write moments to file



numMoments = length(moments{1}.Names);
numOrderOneMoments = sum(moments{1}.Order==1);
numOrderTwoMoments = sum(moments{1}.Order==2);
numOrderThreeMoments = sum(moments{1}.Order==3);
numOrderFourMoments = sum(moments{1}.Order==4);

fName = 'ODEInputMomentsLN';
fid = fopen([fName '.m'], 'w');

fprintf(fid, 'function [dM, Likelihood] = %s(t, M, model)\n\n', fName);
fprintf(fid, 'dM = M;\n');
fprintf(fid, 'update=model.Update;\n');

for i=1:length(model.c)
    if (i==model.InputRateIdx)
        continue;
    end
    fprintf(fid, 'c%d = model.c(%d);\n', i, i);
end



counter = 1;
for k=1:numMoments
    fprintf(fid, '%s = M(%d);\n', moments{1}.Names{k}, counter);
    counter = counter + 1;
end


fprintf(fid, 'Z = GetPieceWiseConstantInput(t, model.TranscriptionRateSequence);');


syms P;
fprintf(fid, '\n');
fprintf(fid, 'if (update==0)\n');

idxCount = 1;


for k=1:numMoments
    
    fprintf(fid, '\tdM(%d) = %s;\n', idxCount, char(moments{1}.dM(k)));
    idxCount = idxCount + 1;
end

fprintf(fid, 'else\n');
fprintf(fid, '\tw = zeros(%d, 1);\n', numSpecies);
fprintf(fid, '\tw(%d) = 1;\n', model.ObservedSpeciesIdx);
fprintf(fid, '\tE=ones(%d, 1);\n', numSpecies);
fprintf(fid, '\tv=model.MeasurementSigma;\n\ty=log(model.Measurement+1);\n');
i = 1;


fprintf(fid, '\n\tMuPrior = zeros(%d, 1);\n', numSpecies);
for k=1:numOrderOneMoments
    fprintf(fid, '\tMuPrior(%d) = %s;\n', k, moments{i}.Names{k});
end

Perms = npermutek(1:numSpecies, 2);
Perms = sort(Perms, 2);
Perms = unique(Perms, 'rows');

startIdx = numOrderOneMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments;
for k=startIdx:endIdx
    fprintf(fid, '\tSPrior(%d, %d) = %s;\n', ...
        Perms(k-numSpecies, 1), Perms(k-numSpecies, 2), moments{i}.Names{k});
    if (numSpecies>1)
        fprintf(fid, '\tSPrior(%d, %d) = %s;\n', ...
            Perms(k-numSpecies, 2), Perms(k-numSpecies, 1), moments{i}.Names{k});
    end
end

startIdx = numOrderOneMoments+numOrderTwoMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments;
for k=startIdx:endIdx
    fprintf(fid, '\tKPrior(%d) = %s;\n', ...
        k-startIdx+1, moments{i}.Names{k});
end

startIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments+numOrderFourMoments;
for k=startIdx:endIdx
    fprintf(fid, '\tUPrior(%d) = %s;\n', ...
        k-startIdx+1, moments{i}.Names{k});
end
%Log[P] + Log[(M + MT + P + S)] - Log[(M + P) (MT + P)]

fprintf(fid, '\tMuPriorLN = 2*log(MuPrior) - 1/2*log(diag(SPrior));\n');
fprintf(fid, '\tSigmaPriorLN = log(SPrior) - log(MuPrior)*E'' - E*log(MuPrior)'';\n');


%fprintf(fid, '\tSigmaPrior = SPrior%d - MuPrior%d*MuPrior%d'';\n', i, i);
fprintf(fid, '\tInvSigmaPriorLN = inv(SigmaPriorLN);\n');
fprintf(fid, '\tInvSigmaLN = w*w'' / v^2 + InvSigmaPriorLN;\n');
fprintf(fid, '\tSigmaPosteriorLN = inv(InvSigmaLN);\n');
fprintf(fid, '\tMuPosteriorLN = SigmaPosteriorLN * (w*y / v^2 + InvSigmaPriorLN * MuPriorLN);\n');
%SigmaPosterior = inv(InvSigma);
%	MuPosterior = SigmaPosterior * (w*y / v^2 + inv(SigmaPrior) * MuPrior);


fprintf(fid, '\tm = w''*MuPriorLN;\n');
fprintf(fid, '\ts = SigmaPriorLN(%d, %d);\n', model.ObservedSpeciesIdx, model.ObservedSpeciesIdx);
%fprintf(fid, '\tLikelihood(%d) =  -1/2*(m - y)^2/(s + v^2) - 1/2*log(s+v^2) ;\n', i);
%fprintf(fid, '\tLikelihood(%d) =  1 - m^2 - m/(2 + v^2) - y + (m*y)/(s + v^2) + y^2/(2*(s + v^2)) - 1/2*log(s + v^2);\n', i);
fprintf(fid, '\tLikelihood = -1/2*(y-m)^2/(v^2 + s) - 1/2*log(1/v^2 + 1/s) - 1/2*log(s) - log(v);\n');
fprintf(fid, '\tif isnan(Likelihood) || isinf(Likelihood)\n');
fprintf(fid, '\t\tLikelihood = -inf;\n');
fprintf(fid, '\t\tMuPosterior = zeros(size(MuPosterior));\n');
fprintf(fid, '\t\tSPosterior = zeros(size(SPosterior));\n');
fprintf(fid, '\t\tSigmaPosterior = zeros(size(SigmaPosterior));\n');
fprintf(fid, '\telse\n');
fprintf(fid, '\tend\n');

%fprintf(fid, '\tPPost(%d) = exp(Likelihood(%d))*P(%d);\n', i, i);

fprintf(fid, '\tMuPosteriorZ=exp(MuPosteriorLN + 1/2*diag(SigmaPosteriorLN));\n');

startIdx = numOrderOneMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments;

for k=startIdx:endIdx
    idxVec = moments{i}.SpeciesIdx{k};
    fprintf(fid, '\tSPosteriorZ(%d, %d) = exp(MuPosteriorLN(%d) + MuPosteriorLN(%d) + ...\n', ...
        idxVec(1), idxVec(2), idxVec(1), idxVec(2));
    fprintf(fid, '\t\t1/2*(SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d)+...\n',...
        idxVec(1), idxVec(1), idxVec(1), idxVec(2));
    fprintf(fid, '\t\tSigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d)));\n',...
        idxVec(2), idxVec(1), idxVec(2), idxVec(2));
    
    fprintf(fid, '\tSPosteriorZ(%d, %d) = SPosteriorZ(%d, %d);\n', idxVec(2), idxVec(1), idxVec(1), idxVec(2));
end

fprintf(fid, '\n\tMuPosterior=MuPosteriorZ;\n');
fprintf(fid, '\tSPosterior=SPosteriorZ;\n\n');

startIdx = numOrderOneMoments+numOrderTwoMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments;
%E^(m1 + m2 + m3 + 1/2 (s11 + s12 + s13 + s21 + s22 + s23 + s31 + s32 + s33))
for k=startIdx:endIdx
    idxVec = moments{i}.SpeciesIdx{k};
    fprintf(fid, '\tKPosterior(%d) = exp(MuPosteriorLN(%d) + MuPosteriorLN(%d) + MuPosteriorLN(%d) + ...\n', ...
        k-startIdx+1, idxVec(1), idxVec(2), idxVec(3));
    fprintf(fid, '\t\t1/2*(SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) +...\n',...
        idxVec(1), idxVec(1), idxVec(1), idxVec(2), idxVec(1), idxVec(3));
    fprintf(fid, '\t\tSigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + ...\n',...
        idxVec(2), idxVec(1), idxVec(2), idxVec(2), idxVec(2), idxVec(3));
    fprintf(fid, '\t\tSigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d)));\n\n',...
        idxVec(3), idxVec(1), idxVec(3), idxVec(2), idxVec(3), idxVec(3));
end


startIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments+numOrderFourMoments;


%formula for 4th non-central moment from Mathematica

for k=startIdx:endIdx
    idxVec = moments{i}.SpeciesIdx{k};
    fprintf(fid, '\tUPosterior(%d) = exp(MuPosteriorLN(%d) + MuPosteriorLN(%d) + MuPosteriorLN(%d) + MuPosteriorLN(%d) + ...\n', ...
        k-startIdx+1, idxVec(1), idxVec(2), idxVec(3), idxVec(4));
    fprintf(fid, '\t\t1/2*(SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + ...\n',...
        idxVec(1), idxVec(1), idxVec(1), idxVec(2), idxVec(1), idxVec(3), idxVec(1), idxVec(4));
    fprintf(fid, '\t\tSigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) +...\n',...
        idxVec(2), idxVec(1), idxVec(2), idxVec(2), idxVec(2), idxVec(3), idxVec(2), idxVec(4));
    fprintf(fid, '\t\tSigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) +...\n',...
        idxVec(3), idxVec(1), idxVec(3), idxVec(2), idxVec(3), idxVec(3), idxVec(3), idxVec(4));
    fprintf(fid, '\t\tSigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d) + SigmaPosteriorLN(%d, %d)));\n\n',...
        idxVec(4), idxVec(1), idxVec(4), idxVec(2), idxVec(4), idxVec(3), idxVec(4), idxVec(4));
end



fprintf(fid, '\t%% Store values in change vector\n');




for k=1:numSpecies
    idx = k;
    fprintf(fid, '\tdM(%d) = MuPosterior(%d) - MuPrior(%d);\n', ...
        idx, k, k);
end

Perms = npermutek(1:numSpecies, 2);
Perms = sort(Perms, 2);
Perms = unique(Perms, 'rows');

startIdx = numOrderOneMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments;
for k=startIdx:endIdx
    idx = (i-1)*numMoments + k;
    fprintf(fid, '\tdM(%d) = SPosterior(%d, %d) - SPrior(%d, %d);\n', ...
        idx, Perms(k-numSpecies, 1), Perms(k-numSpecies, 2), Perms(k-numSpecies, 1), Perms(k-numSpecies, 2));
end

startIdx = numOrderOneMoments+numOrderTwoMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments;
for k=startIdx:endIdx
    idx = (i-1)*numMoments + k;
    fprintf(fid, '\tdM(%d) = KPosterior(%d)- KPrior(%d);\n', ...
        idx, k-startIdx+1, k-startIdx+1);
end

startIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments+1;
endIdx = numOrderOneMoments+numOrderTwoMoments+numOrderThreeMoments+numOrderFourMoments;
for k=startIdx:endIdx
    idx = (i-1)*numMoments + k;
    fprintf(fid, '\tdM(%d) = UPosterior(%d) - UPrior(%d);\n', ...
        idx, k-startIdx+1, k-startIdx+1);
end


fprintf(fid, 'end\n');
fclose(fid);

odeFun = str2func(fName);


%% create default initial conditions

M0 = zeros(numMoments, 1);

for k=1:numMoments
    
    specIdx = moments{1}.SpeciesIdx{k};
    M0(k) = prod(X0(specIdx));
    
end

infos.DefaultInitialConditions = M0;
infos.MomentSystem = moments;


