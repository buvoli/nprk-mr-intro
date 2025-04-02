function [A, b, name] = NPRK_IMEX_121()
%NP_IMEX_121 first-order IMEX-NPRK Euler method
%   121 => first order, two stage, one-solve

name = 'IMEX-NPRK-1[21]';

A = zeros(2,2,2);
A(2,2,1) = 1;

b = [];

end

