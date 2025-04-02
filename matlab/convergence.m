function [error, leg_str, time] = convergence(methods, problem, Nts, y_ref)
%GETERROR Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 4)
    y_ref = problem.referenceSolution();
end

num_nts = length(Nts);
num_methods = length(methods);

error = zeros(num_nts, num_methods);
time  = zeros(num_nts, num_methods);

for i = 1 : num_nts
    for j = 1 : num_methods
        start_time = tic;
        [y_approx] = methods{j}.solve(problem, Nts(i));
        time(i,j) = toc(start_time);
        error(i,j) = norm(y_approx - y_ref, inf);
    end
end

leg_str = cellfun(@(m) m.name, methods, 'UniformOutput', false);

end
