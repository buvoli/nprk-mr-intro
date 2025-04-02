function [A, b, name] = NP_IMEX_243b_SiSa()
%%NP_IMEX_243b_SiSa 2nd-order, L-Stable, stiffly-accurate SDIRK with 3 implicit solves. 
%  243b_SiSa => second-order, four stage, three implicit solve, Varient b,
%               singly-implicit, stiffly-accurate

name = 'IMEX-NPRK2[43]b-SiSa';

gamma = 0.325754;
fg    = sqrt(1 - 4 * gamma^2 * (gamma * (3 * gamma - 8) + 3));

a321 = (1 - 2 * gamma^2 - fg) / (4 * gamma);
a421 = (-1 + 4 * gamma - 2 * gamma^2 - fg) / (4 * gamma);
a432 = (1 - 2 * gamma^2 + fg) / (4 * gamma);

A = zeros(4,4,4);
A(2,2,1) = gamma;
A(3,2,1) = a321;
A(3,3,2) = gamma;
A(4,2,1) = a421;
A(4,3,2) = a432;
A(4,4,3) = gamma;

b = [];

end

