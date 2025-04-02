function [A, b, c] = ssp33()
%SSP33 Optimal 3rd-order 3 stage SSP method

A = [
    0,      0,      0; 
    1,      0,      0; 
    1/4,    1/4,    0
];
b = [1/6,   1/6,    2/3];
c = sum(A,2);

end