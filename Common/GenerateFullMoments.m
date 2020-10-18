function [odeFun, infos] = GenerateFullMoments(modelIn)


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
        H = subs(model.H, switchRateStr, [switchRateStr '*u']);
        
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
moments{1}.dM = SubstituteMoment(XCond, moments{1}.dM, maxOrder+1);



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
    fprintf('%d\n', j);
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

fName = 'ODEFullMoments';
fid = fopen([fName '.m'], 'w');

fprintf(fid, 'function [dM] = %s(t, M, model)\n\n', fName);
fprintf(fid, 'dM = M;\n');
fprintf(fid, 'update=model.Update;\n');

for i=1:length(model.c)
    %if (i==model.InputRateIdx)
    %    continue;
    %end
    fprintf(fid, 'c%d = model.c(%d);\n', i, i);
end



counter = 1;
for k=1:numMoments
    fprintf(fid, '%s = M(%d);\n', moments{1}.Names{k}, counter);
    counter = counter + 1;
end


fprintf(fid, 'u = GetPieceWiseConstantInput(t, model.InputParams);');


syms P;
fprintf(fid, '\n');
fprintf(fid, 'if (update==0)\n');

idxCount = 1;


for k=1:numMoments
    
    fprintf(fid, '\tdM(%d) = %s;\n', idxCount, char(moments{1}.dM(k)));
    idxCount = idxCount + 1;
end

fprintf(fid, 'else\n');


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


