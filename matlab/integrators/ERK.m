classdef ERK < Integrator
    %ERK Generic explicit RK Integrator
    %   Detailed explanation goes here
    
    properties
        name = ''
    end
    
    properties (SetAccess = protected)
        A;
        b;
        AT;
        bT;
    end
    
    methods
        
        function this = ERK(coeffHandle, name_override)
            %ERK Construct an instance of this class
            %   Detailed explanation goes here
            
            [this.A, this.b] = coeffHandle();
            
            this.AT = transpose(this.A);
            this.bT = this.b(:);

            if(nargin < 2)
                s = functions(coeffHandle);
                this.name = s.function;
            else
                this.name = name_override;
            end 

        end
        
        function yn = step(this, yn, h, problem)
            %STEP advances method by one timestep h
            %   Detailed explanation goes here

            s   = size(this.A, 1);
            dim = length(yn);

            K = zeros(dim, s);            
            for j = 1 : s
                Y_j = yn + K(:,1:j-1) * this.AT(1:j-1, j);
                K(:,j) = h * problem.F(Y_j, Y_j);
            end            
            yn = yn + K * this.bT;
            
        end
                
    end
    
end
