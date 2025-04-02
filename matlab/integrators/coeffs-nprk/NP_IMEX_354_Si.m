function [A, b, name] = NP_IMEX_354_Si()
%%NP_IMEX_354_Si 3rd-order, L-Stable, stiffly accurate, DIRK with 4 implicit stages. 
% 354_Si => third-order, five stage, four implicit solves, singly-implicit

name  = 'IMEX-NPRK3[54]-Si';
gamma = 0.54;

A = zeros(5, 5, 5);

A(2,2,1) = gamma;

A(3,2,1) = 0.1040208587459659;
A(3,3,2) = gamma;

A(4,2,1) = -1.240968174302810;
A(4,3,2) = 	0.4238348297973843;
A(4,4,3) = gamma;

A(5,2,1) = 0.4290344770836952;
A(5,3,2) = -1.082995008615554;
A(5,4,3) = 0.2465116558063914;
A(5,5,4) = gamma;

b = zeros(5, 5);

b(2,1) = -0.3205828811598456;
b(3,2) = 1.009514097875651;
b(4,3) = 0.04458528147075302;
b(5,4) = 0.266483501813441;

end
