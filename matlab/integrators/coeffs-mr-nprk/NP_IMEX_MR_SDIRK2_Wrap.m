function [A, b, name] = NP_IMEX_MR_SDIRK2_Wrap(methodHandle, nu, gamma_selector) %
%NP_IMEX_MR_SDIRK2_WRAP second-order implicitly wrapped MR-NPRK
%   methodHandle - function handle [A,b,c] = methodHandle() that returns
%                  explict method to be wrapped.
%   nu           - final implicit solve Y_s = b + F(Y_s, Y_nu) 
%                  for nu \in {1, ..., shat + 1}.
%  gamma_selector - select implicit coefficient
%                       gamma_selector =  1 -> gamma = (2 + sqrt(2)) / 2
%                       gamma_selector = -1 -> gamma = (2 - sqrt(2)) / 2

[A_hat, b_hat] = methodHandle();
shat = size(A_hat, 1);

if(nargin < 2)
    nu = shat + 1;
end
if(nargin < 3)
    gamma_selector = 1;
end

if(gamma_selector == 1)
    gamma = (2 + sqrt(2)) / 2;
else
    gamma = (2 - sqrt(2)) / 2;
end

name = sprintf('IMEX-MR-NPRK2-SDIRK2-Wrap[%i]', shat);

s = shat + 2;
A = zeros(s,s,s);

% first implicit solve
A(2,2,1) = gamma;
% internal stages
for i = 3 : shat + 1
    A(i,2,1) = A_hat(i - 1, 1);
    for j = 3 : i - 1
        A(i,2,j) = A_hat(i - 1, j - 1);
    end
end
% final implicit solve
A(s, 2, [1 3:shat+1]) = b_hat;
A(s, 2, nu) = b_hat(nu - 1) - gamma;
A(s, s, nu) = gamma;

b = [];

end
