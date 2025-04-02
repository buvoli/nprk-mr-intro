function [A, b, name] = NP_IMIM_232b(transpose)
%NP_IMIM_232b second-order, IMIM-NPRK Midpoint/Crank-Nicolson method
%  232b => second-order, three stage, two implicit solves, variant b

if(nargin < 1)
    transpose = false;
end

name = 'IMIM-NPRK2[32]b';

% IM-IM midpoint / Crank-Nicolson method
A = zeros(3,3,3);
A(2, 2, 1) = 1/2;
A(3, 2, 1) = 1/2;
A(3, 2, 3) = 1/2;

b = [];

if(transpose)
    A = permute(A, [1 3 2]);
end

end

