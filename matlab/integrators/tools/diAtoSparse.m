function [Al, Ad] = diAtoSparse(A, b)
% DIATOSPARSE Converts the rank 3 tensor a_{ijk} and matrix b_{ij} for 
% a diagonally-implicit NP-RK method into a sparse format that describes
% lower triangular (Al) and diagonal entries (Ad).
%
% -- Inputs ---------------------------------------------------------------
%
%   A - rank 3 coefficient tensor a_{ijk} of dimension s x s x s 
%   b - matrix b_{ij} or empty if method is 
%
% -- Outputs --------------------------------------------------------------
%
%   Al - (s x 3) cell array describing lower triangular entries (explicit)
%        Al{i} contains a 3xm array such that:
%
%           Y_i = y_n + h \sum_{j=1}^m Al{i}(3,j) F(Y_{Al{i}(1,j)}, Y_{Al{i}(2,j}) + IMPLICIT-TERMS
%
%   Ad - (s x 3) cell array describing diagonal entries (implicit)
%
%           Y_i = EXPLICIT-TERMS + h \sum_{j=1}^m Ad{i}(j,3) F(Y_{Ad{i}(j,1)}, Y_{Ad{i}(j,2})
%
%   Important:
%     - Ad{i} = [] if Y_i is an explicit stage
%     - Ad{i}(j,1) = i if Y_i is implicit in the first argument of F(y,y)
%     - Ad{i}(j,2) = i if Y_i is implicit in the second argument of F(y,y)
% -------------------------------------------------------------------------

% If method is not stiffly accurate, then add row to A
if(~isempty(b))
    s  = size(A,1);
    Ah = zeros(s+1,s+1,s+1);
    Ah(1:s, 1:s, 1:s) = A;
    Ah(s+1, 1:s, 1:s) = b;
    A = Ah;
end

% -- Convert A to sparse format -------------------------------
s = size(A, 1);
Al = cell(s,1);
Ad = cell(s,1);

for i = 1 : s
    
    % ----> lower diagonals (explicit entries)
    [j, k, val] = find(squeeze(A(i, 1:i-1, 1:i-1))); % extract non zero indices and value
    if(~isempty(j))
        Al{i} = [j, k, val];
    end
    
    % ----> diagonals (implicit entries)
    if(~isempty(find(A(i, i, 1:i-1), 1))) % implicit in arg 1
        
        [~, k, val] = find(A(i, i, 1:i-1));
        j = i * ones(size(k));
        Ad{i} = [j, k, val];
        
    elseif(~isempty(find(A(i, 1:i-1, i), 1))) % implicit in arg 2
        
        [~, k, val] = find(A(i, 1:i-1, i));
        j = i * ones(size(k));
        Ad{i} = [k, j, val];
        
    end
    
end

end


