function [A, b, name] = NP_IMEX_MR_SDIRK3_v1_Wrap(methodHandle, om) %
%NP_IMEX_MR_SDIRK3_WRAP third-order IMEX explicitly wrapped method (Variant
%1)
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

for i = 4 : s2 + 2
    A(i, 2, [1 4:i-1]) = A2(i - 2, 1:(i-3));
end

% special case A_{s+2,3,s2+2} to satisfy b^{2}_i A^{1}_{ij} c_j
A(s2 + 2, 3, s2 + 1) = (1 - 3 * c1(1)) / ( 6 * b2(s2) * (c1(2) - c1(1)));
A(s2 + 2, 2, s2 + 1) = A2(s2, s2 - 1) - A(s2 + 2, 3, s2 + 1);

% -- output stage ---------------------------------------------------------

A(s, 2, [1 4:s2+2]) = b2;

% special case a_{s,3,1} and a_{s,3,s-1}
eta            = outputCoeffF2(A1, b1, c1, A2, b2, c2);
A(s, 3, 1)     = eta;
A(s, 3, s - 1) = A1(3, 2) - eta;

A(s, 2, 1)     = b2(1) - A(s, 3, 1);
A(s, 2, s - 1) = b2(end) - A(s, 3, s-1) - A1(3, 3);
A(s, s, s - 1) = A1(3, 3);

b = []; % empty b for stiffly accurate method

end

function eta = outputCoeffF2(A1, b1, c1, A2, b2, c2)
% order condition for c_i A_{ij} c_j = 1 / 3 for G(w) = F(Y2,w) in output

    n = 1 / 3 - A1(3,3) * c2(end) + c1(1) * (A1(3,3) * c2(end) - 1 / 2) - A1(3,2) * (c1(2) - c1(1)) * c2(end);
    d = c2(end) * (c1(1) - c1(2));
    eta = n / d;

end
