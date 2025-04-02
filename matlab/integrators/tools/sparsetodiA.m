function [A, b] = sparsetodiA(Al, Ad)
%SPARSETODIA Convert sparse diagonally implicit NP-RK format to the rank 3 tensor a_{ijk} and matrix b_{ij}
%   Detailed explanation goes here

s = length(Al);
A = zeros(s,s,s);

for i = 1 : s    
    for j = 1 : size(Al{i},1)
        A(i, Al{i}(j, 1), Al{i}(j, 2)) = Al{i}(j, 3);
    end
    for j = 1 : size(Ad{i},1)
        A(i, Ad{i}(j, 1), Ad{i}(j, 2)) = Ad{i}(j, 3);
    end
end

if(isempty(Ad{s})) % non-stiffly accurate method
    b = squeeze(A(s, 1 : s - 1, 1 : s - 1));
    A = A(1 : s - 1, 1 : s - 1, 1 : s - 1);
else % stiffly accurate method
    b = [];
end

end
