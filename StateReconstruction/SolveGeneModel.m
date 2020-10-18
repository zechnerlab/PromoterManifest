function output = SolveGeneModel(t, model, config, transcrSequence, M0)


if (nargin<5)
M0 = model.InputODEInfos.infos.DefaultInitialConditions;
end
transT = transcrSequence.t;
Transcr = transcrSequence.Value;

tStart = t(1);
tEnd = t(end);


bIdx = find(transT>tStart);
lIdx = find(transT<tEnd);
loopIdx = intersect(bIdx, lIdx);

tVec = [tStart, transT(loopIdx), tEnd];
transcrVal = [Transcr(bIdx(1)-1), Transcr(loopIdx)];


currT = tVec(1);

MTot = M0;
tTot = currT;

opts = {'RelTol', 0.0001, 'AbsTol', 0.0001};

output.valid = 1;

for j=2:length(tVec)
    
    nextT = tVec(j);

    if (nextT-currT)<eps
       continue; 
    end
    
    numSteps = max(3, ceil((nextT - currT) / config.ODEStepSize));
    tGrid = linspace(currT, nextT, numSteps);
    
    try
        model.Z = transcrVal(j-1);
        model.Update = 0;
        if (config.ODEType==1)
            odeFun = @(t, y) model.InputODEInfos.odeHandle(t, y, model);
            M = ode2(odeFun, tGrid, M0);
            
        elseif (config.ODEType==2)
            [~, M] = ode45(model.InputODEInfos.odeHandle, tGrid, M0, opts, model);
        end
        
        M = M';
        
    catch msg
        fprintf('something went wrong!\n');
        output.valid = 0;
        return;
    end
    
    MTot = [MTot, M(:, 2:end)];
    tTot = [tTot, tGrid(2:end)];
    
    M0 = M(:, end);
    
    currT = nextT;
    
end

output.M = MTot;
output.t = tTot;


