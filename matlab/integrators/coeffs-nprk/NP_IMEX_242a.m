function [A, b, name] = NP_IMEX_242a()
%NP_IMEX_242a second-order, L-stable IMEX-NPRK (global stability optimum)
%  232 => second-order, four stage, two implicit solve

name  = 'IMEX-NPRK2[42]a';
gamma = 1 + 1 / sqrt(2);

A = zeros(4, 4, 4);

A(2,2,1) = gamma;
A(3,2,1) = (1 / 42) * (26 - 3 * sqrt(2));
A(4,2,1) = (1 / 42) * (-20 - 23 * sqrt(2));
A(4,4,3) = gamma;

b = zeros(4,4);

b(2,1) = (1 / 94) * (16 - 9 * sqrt(2));
b(4,3) = (3 / 94) * (26 + 3 * sqrt(2));

end

