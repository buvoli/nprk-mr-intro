classdef NPRK < Integrator
    %NPRK A Generic Diagonally-Implicit NPRK Integrator
    %
    %       Y_i = y_n + \sum_{j=1}^{i} \sum_{k=1}^{i} a_{ijk} F(Y_j,Y_k)
    %       y_{n+1} = y_n + \sum_{j=1}^{i}
    %
    %   with the additional condition that each stage is implicit in only
    %   one argument:
    %
    %       a_{i,i,k} \ne 0 \iff a_{i,j,i} = 0
    %       a_{i,j,i} \ne 0 \iff a_{i,i,k} = 0

    properties
        name = ''
    end
    
    properties (SetAccess = protected)
        Al;
        Ad;
    end
    
    methods
        
        function this = NPRK(coeffHandle)
            %SemiImplicitRK2B Construct an instance of this class
            %   Detailed explanation goes here
            
            [this.Al, this.Ad, this.name] = convertTableau(coeffHandle, 'sparse');
            
        end
        
        function yn = step(this, yn, h, problem)
            %STEP advances method by one timestep h
            %   Detailed explanation goes here
            
            s = length(this.Al);
            dim = length(yn);
            
            Y = zeros(dim, s); % store stage values Y_i
            hF = cell(s,s);     % stores stage derivative F(Y_j,Y_k)
            
            for i = 1 : s
                
                Al_i = this.Al{i};
                Ad_i = this.Ad{i};
                
                % ----> form explicit vector 
                b = yn;
                
                for l = 1 : size(Al_i,1)
                    
                    j = Al_i(l, 1);
                    k = Al_i(l, 2);
                    a_jk = Al_i(l, 3);
                    
                    if( isempty(hF{j,k}) ) % evaluate and store F(Y_j,Y_k)
                        hF{j,k} = h * problem.F(Y(:,j), Y(:,k));
                    end
                    
                    b = b + a_jk * hF{j, k};
                    
                end
                
                % ----> implicit solve
                
                if(isempty(Ad_i)) % explicit method
                    
                    Y(:,i) = b;
                    
                else % implicit method
                    
                    if(all(Ad_i(:,1) == i)) % implicit in first arg
                    
                        Y(:,i) = problem.invertFComp(b, h * Ad_i(:,3), Y(:, Ad_i(:,2)), yn, 1);
                        
                    else % implicit in second first arg
                    
                        Y(:,i) = problem.invertFComp(b, h * Ad_i(:,3), Y(:, Ad_i(:,1)), yn, 2);
                        
                    end
                    
                end
                
            end
            
            yn = Y(:,end);
            
        end
                
    end
    
end
