function [A, b, name] = NP_IMEX_354_Sa()
%%NP_IMEX_354_Sa 3rd-order, L-Stable, SDIRK with 4 implicit stages (gamma = 54/100). 
%  354_Sa => third-order, five stage, four implicit solves, stiffly-accurate

name = 'IMEX-NPRK3[54]-Sa';

A = zeros(5, 5, 5);
                
A(2,2,1) =   1;

A(3,2,1) = - 2 / 3;
A(3,3,2) =   2 / 3;

A(4,2,1) =   5 / 12;
A(4,3,2) = - 5 / 12;
A(4,4,3) =   1 / 2;

A(5,2,1) = - 1 / 2;
A(5,3,2) =   1 / 6;
A(5,4,3) =   2 / 3;
A(5,5,4) =   2 / 3;

b = [];

end
