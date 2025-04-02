function [A, b, name] = NP_IMEX_232a()
%NP_IMEX_232a second-order, L-stable IMEX-NPRK (global stability optimum)
%  232 => second-order, three stage, two implicit solve 

name  = 'IMEX-NPRK2[32]a';
gamma = (1 + 1/sqrt(2));

A = zeros(3, 3, 3);

A(2,2,1) = gamma;
A(3,2,1) = (-2 - 3/sqrt(2));
A(3,3,2) = gamma;

b = zeros(3,3);

b(2,1) = 1 / sqrt(2);
b(3,2) = 1 - 1 / sqrt(2);

end
