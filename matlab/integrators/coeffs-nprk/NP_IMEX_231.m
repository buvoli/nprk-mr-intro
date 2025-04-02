function [A, b, name] = NP_IMEX_231()
%NP_IMEX_231 second-order, IMEX-NPRK midpoint method
%  231 => second-order, three stage, one implicit solve

name = 'IMEX-NPRK2[31]';

A = zeros(2,2,2);

A(2,2,1) = 1/2;

b = zeros(2,2);

b(2,2) = 1;

end

