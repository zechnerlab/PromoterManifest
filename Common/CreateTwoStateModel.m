function model = CreateTwoStateModel(on, off, leak, transcr)

    model.NumStates = 2;

    Q = zeros(model.NumStates, model.NumStates);
    Q(1, 2) = off;
    Q(2, 1) = on;
    Q = Q - diag(sum(Q, 1));
    Z = [leak; transcr];
    
    model.NumStates = 2;
    model.Z = Z;
    model.Q = Q;
    
    model.P0 = [1;0]
end