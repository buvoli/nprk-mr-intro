function [nA, nb, nc] = nStepTableau(arg1, arg2, arg3)
%NSTEPTABLEAU returns a tablue cooresponding to n steps of method (A,b,c)
% [nA, nB] = nStepTablea(methodHandle, n)
%       methodHandle - returns [A,b,c] of RK method
%       n - number of steps
% [nA, nb] = nStepTablea(A, b, n)
%       A - Butcher matrix
%       b - weight vector
%       n - number of steps

if(nargin == 2)
    methodHandle = arg1; n = arg2;
    [A, b] = methodHandle();
else
    A = arg1; b = arg2; n = arg3;
end

s  = size(A, 1);
M1 = eye(n);
M2 = tril(ones(n), -1);
B  = repmat(b(:).', s, 1);

nA = (kron(M1, A) + kron(M2, B)) / n;
nb = repmat(b(:), n, 1) / n;
nc = sum(nA, 2);

end