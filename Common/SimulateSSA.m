% Wrapper function for the SimulateMarginalProcess function to simulate a 
% standard SSA path.
function [X, t, r, G] = SimulateSSA(X0, c, Pre, Post, N, startTime, T)


    [X, t, r, G] = ...
        SimulateMarginalProcess(X0, c, Pre,...
        Post, N, T, 1, [0, T], c(1)*[1, 1, 1], [T/2, T], startTime, ...
        [], [], [], [], []);

end