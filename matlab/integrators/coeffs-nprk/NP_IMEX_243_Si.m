function [A, b, name] = NP_IMEX_243_Si()
%%NP_IMEX_243_Si 2nd-order, L-Stable, stiffly-accurate SDIRK with 3 implicit solves. 
%  243_Si => second-order, four stage, three implicit solve, singly-implicit

name = 'IMEX-NPRK2[43]-Si';

gamma = 0.553658;
b32 = -0.0054849;
b43 = 0.237378;

a321 = (1 - 2 * gamma * (b32 + b43)) / (2 * b43);
a432 = (gamma * (-1 - 2 * (-2 + gamma) * gamma)) / (-1 + 2 * gamma * (b32 + b43)); 
b21  = 1 - b32 - b43;
a421 = (1 / 2) * ((b32 * (-1 + 2 * gamma * b32)) / b43^2 + (1 + 2 * gamma * (-1 + b32)) / b43 + (2 * gamma * (1 + 2 * (-2 + gamma) * gamma)) / (-1 + 2 * gamma * (b32 + b43)));

% -------------------------------------------------------------------------

A = zeros(4,4,4);
A(2,2,1) = gamma;
A(3,2,1) = a321;
A(3,3,2) = gamma;
A(4,2,1) = a421;
A(4,3,2) = a432;
A(4,4,3) = gamma;

b = zeros(4,4);
b(2,1) = b21;
b(3,2) = b32;
b(4,3) = b43;

end

