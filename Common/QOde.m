function [dP] = QOde(t, P, Q, Z, update)

dP = zeros(size(P));

M = Z'*P(1:end-1);


dP(1:end-1) = (Q - diag(Z - M))*P(1:end-1);
dP(end) = M;

end