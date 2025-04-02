classdef Problem < handle
    %PROBLEM Abstract problem class
    
    properties (Abstract)
        tspan
        initial_condition
    end
    
    methods (Abstract)
        F(this, y1, y2) % evaluates full right hand side F(y,y)
        invertFComp(this, x, a, y, U0, component) 
        % If component = 1, returns solution U to U = x + \sum_{j=1}^n a_j H(U, y(j)). 
        % If component = 2, returns solution U to U = x + \sum_{j=1}^n a_j H(y(j), U).
        % U0 is initial guess.
    end
    
    methods 
    
        function y_ref = referenceSolution(this)
            options = odeset('RelTol', 5e-14, 'AbsTol', 5e-14);
            [~,ys_ref] = ode45(@(t,y) this.F(y,y),[this.tspan(1) mean(this.tspan) this.tspan(end)],this.initial_condition,options); 
            y_ref = ys_ref(end, :).';
        end

        function [ts, ys_ref] = referenceSnapshots(this, num_snapshots)
            
            if(nargin < 2)
                num_snapshots = 500;
            end
               
            options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
            ode_tspan = linspace(this.tspan(1), this.tspan(end), num_snapshots);
            [ts, ys_ref] = ode45(@(t,y) this.F(y,y), ode_tspan, this.initial_condition, options); 

            if(nargout == 1)
                ts = ys_ref;
            end

        end
    
    end
    
end