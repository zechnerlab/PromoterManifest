clear all;
close all;

c1 = 0.0011;
c2 = 100;
c3 = 1;

Q = [0  c1 0 
     c1 0  c2
     0  c2 0];

Z = [0, 0.1, 0.1]; 
 
Q = Q - diag(sum(Q));

[u, v, w] = eig(Q);

%[u] = null(Q);

%[u, v] =svd(Q');

U = u;
iU = inv(U);

v(1) = -inf;