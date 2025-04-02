classdef Integrator < handle
    %INTEGRATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract)
        name
    end
    
    methods
        
        function [y] = solve(this, problem, num_steps)
            %INTEGRATOR Construct an instance of this class
            %   Detailed explanation goes here
            
            h = diff(problem.tspan([1, end])) / num_steps;
            y = problem.initial_condition;
            
            for i = 1 : num_steps
                y = this.step(y, h, problem);

                % --> emergency exit conditions
                if(any(isnan(y(:))) || any(isinf(y(:))))
                    y = NaN;
                    warning('Integrator: Emergency exit at step %i', i);
                    break;
                end

            end
             
        end
        
        
        function [ts, ys] = solveAndStore(this, problem, num_steps)
            %INTEGRATOR Construct an instance of this class
            %   Detailed explanation goes here
            
            h = diff(problem.tspan([1, end])) / num_steps;
            y0 = problem.initial_condition;
            
            ts = linspace(problem.tspan(1), problem.tspan(end), num_steps + 1);
            ys = zeros(length(y0), num_steps + 1);
            ys(:,1) = y0;
            
            for i = 1 : num_steps
                ys(:,i+1) = this.step(ys(:,i), h, problem);
            end

            if(nargout == 1)
                ts = ys;
            end
             
        end
        
    end
    
    methods (Abstract)
        step(this, yn, h, problem);
    end
    
end

