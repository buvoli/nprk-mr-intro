function [A, b, name] = NP_IMIM_132(transpose)
%NP_IMIM_132 first-order, IMIM-NPRK method
%  132a => first-order, three stage, two implicit solves

if(nargin < 1)
    transpose = false;
end

name = 'IMIM-NPRK1[32]';

A = zeros(3,3,3);
A(2,2,1) = 1;
A(3,2,3) = 1;

b = [];

if(transpose)
    A = permute(A, [1 3 2]);
end

end

