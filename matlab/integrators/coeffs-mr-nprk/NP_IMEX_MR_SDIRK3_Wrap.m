function [A, b, name] = NP_IMEX_MR_SDIRK3_Wrap(methodHandle, om, nu, es_type, os_type) %
%NP_IMEX_MR_SDIRK3_WRAP third-order IMEX explicitly wrapped method
%   methodHandle - function handle [A,b,c] = methodHandle() that returns
%                  explict method to be wrapped.
%   om           - second implicit solve Y_3 = b + F(Y_3, Y_om)
%                  for om \in {1,2}.
%   nu           - final implicit solve Y_s = b + F(Y_s, Y_nu) 
%                  for nu \in {1, ..., shat + 2}.
%   es_type      - Explicit stage sparsity. Options are "f2" or "f3" or "f23"
%                   f2: function evaluations in explicit stages involve
%                       G(w) = F(Y2, w) 
%                   f3: function evaluations in explicit stages involve
%                       G(w) = F(Y3, w) 
%                   f23: function evaluations in explicit stages involve
%                       G(w) = (1-δ) F(Y2, w) + δ F(Y3, w)
%   os_type      - Output stage sparsity. Options are "f2" or "f3" or "f23"
%                   f2: function evaluations in final implicit stage involve
%                       G(w) = F(Y2, w) 
%                   f3: function evaluations in final implicit stage involve
%                       G(w) = F(Y3, w) 
%                   f23: function evaluations in final implicit stage involve
%                       G(w) = (1-δ) F(Y2, w) + δ F(Y3, w)
%
% REMARK: code is written to allow for further configuration of explicit
% stages through parameters u,w and the output stage through parameters
% alpha, beta. The values of these are hard-coded to avoid more input
% arguments, however the code is writen generically and they can be changed
% manually, or added as additinal arguments

[A1, b1, c1] = SDIRK3_Sa();     % underlying implicit method (reduced)
[A2, b2, c2] = methodHandle();  % underlying explicit method (reduced)
s2           = size(A2, 1);
name         = sprintf('IMEX-MR-NPRK3-SDIRK3-Wrap[%i]', s2);

if(nargin < 2)
    om = 2; % second implicit solve Y_3 = b + F(Y_3, Y_om)
end

if(nargin < 3)
    nu = s2 + 2; % final implicit solve Y_s = b + F(Y_s, Y_nu)
end

if(nargin < 4)
    es_type = 'f23'; % 'f2', 'f3', 'f23'
end

if(nargin < 5)
    os_type = 'f23';  % 'f2', 'f3', 'f23'
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
A(2,2,1)  = A1(1,1);
A(3,2,1)  = A1(2,1) - a322;
A(3,2,2)  = a322;
A(3,3,om) = A1(2,2);

% -- explicit stages ------------------------------------------------------
switch(es_type)
    
    case 'f2' % all but one explicit stage only involves evaluation of G(w) = F(Y2,w)

        for i = 4 : s2 + 2
            A(i, 2, 1) = A2(i - 2, 1);
            for j = 4 : i - 1
                A(i, 2, j) = A2(i - 2, j - 2);
            end
        end

        % special case A_{u,3,w} to satisfy b^{2}_i A^{1}_{ij} c_j
        u = s2 + 2;
        w = s2 + 1;
        A(u, 3, w) = (1 - 3 * c1(1)) / ( 6 * b2(u - 2) * (c1(2) - c1(1)));
        A(u, 2, w) = A2(u - 2, max(1, w - 2)) - A(u, 3, w);
        
    case 'f3' % all but one explicit stage only involves evaluation of G(w) = F(Y3,w)
    
        for i = 4 : s2 + 2
            A(i, 3, 1) = A2(i - 2, 1);
            for j = 4 : i - 1
                A(i, 3, j) = A2(i - 2, j - 2);
            end
        end

        % special case A_{u,3,w} to satisfy b^{2}_i A^{1}_{ij} c_j
        u = s2 + 2;
        w = s2 + 1;
        A(u, 3, w) = (1 - 3 * c1(2)) / ( 6 * b2(u - 2) * (c1(2) - c1(1))) + A2(u-2, max(1, w - 2));
        A(u, 2, w) = A2(u - 2, max(1, w - 2)) - A(u, 3, w);
    
    case 'f23' % all explicit stage evaluate of G(w) = \delta F(Y2,w) + (1-\delta) F(Y3,w)

        delta = - 0.36350683689006809456061097715969; % = (1 - 3 * c1(1))/(3 * (c1(2) - c1(1)));
        
        for i = 4 : s2 + 2
            A(i, 2, 1) = (1 - delta) * A2(i - 2, 1);
            A(i, 3, 1) = (delta) * A2(i - 2, 1);
            for j = 4 : i - 1
                A(i, 2, j) = (1 - delta) * A2(i - 2, j - 2);
                A(i, 3, j) = (delta) * A2(i - 2, j - 2);
            end
        end

end

% -- output stage ---------------------------------------------------------

A(s, s, nu) = A1(3,3);

switch(os_type)
     
    case 'f2'
        
         A(s, 2, [1 4:s2+2]) = b2;
         A(s, 2, nu) = A(s, 2, nu) - A1(3,3); % adjust implicit stage
         
         % special case a_{s,3,alpha} and a_{s,3,beta}
         % (enforce M1 and order condition c_i A_{ij} c_j = 1/3)
         alpha = s2 + 1; 
         beta  = s2 + 2;

         alpha = 1; 
         beta  = s2 + 2;

         A(s, 3, alpha) = outputCoeffF2(A1, c1, A2, c2, alpha, beta, nu);
         A(s, 3, beta)  = A1(3,2) - A(s, 3, alpha);
         
         A(s, 2, alpha) = A(s, 2, alpha) - A(s, 3, alpha);
         A(s, 2, beta)  = A(s, 2, beta) - A(s, 3, beta);

     case 'f3'

         A(s, 3, [1 4:s2+2]) = b2;
         A(s, 2, nu) = A(s, 2, nu) - A1(3,3); % adjust implicit stage
         
         % special case a_{s,3,alpha} and a_{s,3,beta}
         % (enforce M1 and order condition c_i A_{ij} c_j = 1/3)
         alpha = 1; 
         beta  = s2 + 2;

         A(s, 3, alpha) = outputCoeffF3(A1, b1, c1, A2, b2, c2, alpha, beta, nu);
         A(s, 3, beta)  = A1(3,2) - A(s, 3, alpha) - (1 - b2(max(1,alpha-2)) - b2(max(1,beta - 2)));
         
         A(s, 2, alpha) = A(s, 2, alpha) + b2(max(1,alpha-2)) - A(s, 3, alpha);
         A(s, 2, beta)  = A(s, 2, beta) + b2(max(1,beta-2)) - A(s, 3, beta);
    
     case 'f23'

         delta = - 0.36350683689006809456061097715969; % = (1 - 3 * c1(1))/(3 * (c1(2) - c1(1)));
         
         A(s, 2, [1 4:s2+2]) = (1 - delta) * b2;
         A(s, 3, [1 4:s2+2]) = delta * b2;
         
         % special case a_{s,3,alpha} and a_{s,3,beta}
         % (enforce M1 and order condition c_i A_{ij} c_j = 1/3)
         alpha = 1; 
         beta  = s2 + 2;

         [eta, X1, X2, X3] = outputCoeffF23(A1, b1, c1, A2, b2, c2, alpha, beta, nu, delta);
         A(s, 3, alpha) = eta;
         A(s, 3, beta)  = -eta + X1;

         A(s, 2, alpha) = b2(max(1,alpha-2)) - A(s, 3, alpha); % = -eta + X2
         A(s, 2, beta)  = b2(max(1,beta-2)) - A(s, 3, beta);   % eta + X3
         A(s, 2, nu)    = A(s, 2, nu) - A1(3,3);               % adjust implicit stage

end

b = []; % empty b for stiffly accurate method

end

function eta = outputCoeffF2(A1, c1, A2, c2, alpha, beta, nu)
% order condition for c_i A_{ij} c_j = 1 / 3 for G(w) = F(Y2,w) in output

    al_i = max(1, alpha - 2);
    be_i = max(1, beta - 2);
    nu_i = max(1, nu - 2);

    n = 1 / 3 - A1(3,3) * c2(nu_i) + c1(1) * (A1(3,3) * c2(nu_i) - 1 / 2) - A1(3,2) * (c1(2) - c1(1)) * c2(be_i);
    d = (c1(2) - c1(1)) * (c2(al_i) - c2(be_i));
    eta = n / d;

end

function eta = outputCoeffF3(A1, b1, c1, A2, b2, c2, alpha, beta, nu)
% order condition for c_i A_{ij} c_j = 1 / 3 for G(w) = F(Y3,w) in output

    al_i = max(1, alpha - 2);
    be_i = max(1, beta - 2);
    nu_i = max(1, nu - 2);

    % Long Formula (equivalent to the one below)
    % t1 = c1(1) * ( b2(al_i) * c2(al_i) + (b2(be_i) - (A1(3,2) - (1 - b2(al_i) - b2(be_i)))) * c2(be_i) - A1(3,3) * c2(nu_i));
    % t2 = c1(2) * ( 1 / 2 - b2(al_i) * c2(al_i) - b2(be_i) * c2(be_i) + (A1(3,2) - (1 - b2(al_i) - b2(be_i))) * c2(be_i));
    % t3 = A1(3,3) * c2(nu_i);
    % 
    % eta = (1/3 - t1 - t2 - t3) / ((c1(2) - c1(1)) * (c2(al_i) - c2(be_i)));
 
    n11  = 1 / 3 - A1(3,3) * c2(nu_i);
    n12  = - c1(1) * ( b2(al_i) * c2(al_i) + b2(be_i) * c2(be_i) - A1(3,3) * c2(nu_i) );
    n13  = - c1(2) * ( 1 / 2 - b2(al_i) * c2(al_i) - b2(be_i) * c2(be_i));
    n1   = n11 + n12 + n13;
    d1   = (c1(2) - c1(1)) * (c2(al_i) - c2(be_i));

    n2 = -1 * (A1(3,2) - (1 - b2(al_i) - b2(be_i))) * c2(be_i);
    d2 = (c2(al_i) - c2(be_i));

    eta = n1 / d1 + n2 / d2;

end

function [eta, X1, X2, X3, X4] = outputCoeffF23(A1, b1, c1, A2, b2, c2, alpha, beta, nu, delta)
% order condition for c_i A_{ij} c_j = 1 / 3 for G(w) = (1 - delta) F(Y2,w) + delta F(Y3,w) in output

    al_i = max(1, alpha - 2);
    be_i = max(1, beta - 2);
    nu_i = max(1, nu - 2);

    X1 = A1(3,2) - delta * (1 - b2(al_i) - b2(be_i));
    X2 = b2(al_i);
    X3 = b2(be_i) - X1;
    X4 = (1/2) - (b2(al_i) * c2(al_i) + b2(be_i) * c2(be_i));
    
    num = 1 / 3 - c2(nu_i) * A1(3,3) - c1(1) * ((1 - delta) * X4 + c2(al_i) * X2 + c2(be_i) * X3 - A1(3,3) * c2(nu_i)) - c1(2) * ( delta * X4 + c2(be_i) * X1 );
    den = (c1(2) - c1(1)) * (c2(al_i) - c2(be_i));
    
    eta = num / den;

end
