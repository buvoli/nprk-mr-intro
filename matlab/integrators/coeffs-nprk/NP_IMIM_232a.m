function [A, b, name] = NP_IMIM_232a(transpose)
%NP_IMIM_232a second-order, IMIM-NPRK midpoint method
%  232a => second-order, three stage, two implicit solves, variant a

if(nargin < 1)
    transpose = false;
end

name = 'IMIM-NPRK2[32]a';

% IM-IM Midpoint Method
A = zeros(3,3,3);
A(2,2,1) = 1/2;
A(3,2,3) = 1/2; 

b = zeros(3,3);
b(2,3) = 1;

if(transpose)
    A = permute(A, [1 3 2]);
    b = permute(b, [2 1]);
end

end

