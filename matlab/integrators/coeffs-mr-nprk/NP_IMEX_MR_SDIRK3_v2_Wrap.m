function [A, b, name] = NP_IMEX_MR_SDIRK3_v2_Wrap(methodHandle, om) %
%NP_IMEX_MR_SDIRK3_WRAP third-order IMEX explicitly wrapped method (Variant
%2)
%   methodHandle - function handle [A,b,c] = methodHandle() that returns
%                  explict method to be wrapped.
%   om           - second implicit solve Y_3 = b + F(Y_3, Y_om)

[A1, b1, c1] = SDIRK3_Sa();     % underlying implicit method (reduced)
[A2, b2, c2] = methodHandle();  % underlying explicit method (reduced)
s2           = size(A2, 1);
name         = sprintf('IMEX-MR-NPRK3-SDIRK3-Wrap[%i]', s2);

if(nargin < 2)
    om = 2; % second implicit solve Y_3 = b + F(Y_3, Y_om)
end

s = s2 + 3;
A = zeros(s, s, s);

if(om == 1)
    a322 = 0.1825367204468751798095174531383;
elseif(om == 2)
    a322 = - 0.253329801061583819606501998055;
else
    error('invalid omega');
end

% -- first two implicit solves --------------------------------------------
A(2, 2, 1)  = A1(1, 1);
A(3, 2, 1)  = A1(2, 1) - a322;
A(3, 2, 2)  = a322;
A(3, 3, om) = A1(2, 2);

% -- explicit stages ------------------------------------------------------

delta = - 0.36350683689006809456061097715969; % = (1 - 3 * c1(1))/(3 * (c1(2) - c1(1)));

for i = 4 : s2 + 2
    A(i, 2, [1 4:i-1]) = (1 - delta) * A2(i - 2, 1:(i-3));
    A(i, 3, [1 4:i-1]) = (delta) * A2(i - 2, 1:(i-3));    
end

% -- output stage ---------------------------------------------------------

A(s, 2, [1 4:s2+2]) = (1 - delta) * b2;
A(s, 3, [1 4:s2+2]) = delta * b2;

% special case a_{s,3,1} and a_{s,3,s-1}
eta            = outputCoeffF23(A1, b1, c1, A2, b2, c2, delta);
A(s, 3, 1)     = eta;
A(s, 3, s - 1) = A1(3, 2) - eta - delta * (1 - b2(1) - b2(end));

A(s, 2, 1)     = b2(1) - A(s, 3, 1);
A(s, 2, s - 1) = b2(end) - A(s, 3, s-1) - A1(3, 3);
A(s, s, s - 1) = A1(3, 3);

b = []; % empty b for stiffly accurate method

end

function [ eta ] = outputCoeffF23(A1, b1, c1, A2, b2, c2, delta)
% order condition for c_i A_{ij} c_j = 1 / 3 for G(w) = (1 - delta) F(Y2,w) + delta F(Y3,w) in output

    X1 = A1(3,2) - delta * (1 - b2(1) - b2(end));
    X3 = b2(end) - X1;
    X4 = (1/2) - (b2(end) * c2(end));
    
    N1 = ((1 - delta) * X4 + c2(end) * X3 - A1(3,3) * c2(end));
    N2 = (delta * X4 + c2(end) * X1);
    
    num = 1 / 3 - c2(end) * A1(3,3) - c1(1) * N1 - c1(2) * N2;
    den = (c1(1) - c1(2)) * (c2(end));
    
    eta = num / den;

end
