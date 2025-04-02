classdef Burgers < Problem
    %BURGERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Generic properties
        tspan = [0 6/10]
        initial_condition
        % Spatial grid
        Nx    = 200         % number of internal gridpoints
        xspan = [-2, 2]     % spatial domain range
        xs                  % spatial grid
        conservative = 0;   % 0 - non conservative uu_x, 1 - conservative (u^2/2)_x
        % PDE parameters
        epsilon = 0;
        arg_lock_enabled = true;
    end
    
    properties (SetAccess = protected)
        % operators
        DXX
        DX
        I
        % argument lock
        locked_arg_count = 0;  % number of stored locked args
        locked_arg_inds  = []; % stores argument index for ith locked arg
        locked_args      = []; % stores locked argument data    
    end
    
    methods
        
        function this = Burgers()
            %BURGERS Construct an instance of this class
            %   Detailed explanation goes here
            this.reset();
        end
        
        function rhs = F(this, y_imp, y_exp)
            %F ODE right-hand-side
            if(this.conservative == 1)
                rhs = this.epsilon * this.DXX * y_imp + this.DX * (y_exp .* y_imp) / 2;
            else
                rhs = this.epsilon * this.DXX * y_imp + y_exp .* (this.DX * y_imp);
            end
        end
        
        function y = invertFComp(this, b, as, ys, y0, component)
            %INVERTIMP invert system involving right-hand-side
            % If component = 1, returns solution U to U = b + \sum_{j=1}^n as(j) H(U, ys(j))
            % If component = 2, returns solution U to U = b + \sum_{j=1}^n as(j) H(ys(j), U)
            % U0 is initial guess.
            
            switch(component)
                case 1
                    y = invertImp1(this, b, as, ys, y0);
                case 2
                    y = invertImp2(this, b, as, ys, y0);
            end
            
        end
        
        function y = invertImp1(this, x, as, ys, y0)
            %INVERTIMP solves the equation y = x + \sum_{j=1}^n a_j H(y,y_j)            
            
            c_diff = this.epsilon * sum(as); 
            
            dim = length(x);
            c_adv = ys * as(:); % vector \sum_{j=1}^n a(j) y(:,j)
            c_adv = spdiags(c_adv, 0, dim, dim);
            
            % conservative nonlinearity (u^2/2)_x -> D_x * diag(u_exp) * u_imp
            if(this.conservative == 1)
                y = (this.I - c_diff * this.DXX - this.DX * c_adv / 2) \ x;
            else % nonconservative nonlinearity u u_x -> diag(u_exp) * D_x * u_imp
                y = (this.I - c_diff * this.DXX - c_adv * this.DX) \ x;
            end
            
        end
        
        function y = invertImp2(this, x, as, ys, y0)
            %INVERTIMP solves the equation y = x + \sum_{j=1}^n a_j H(y,y_j)            
            
            dim = length(x);
            c_diff = this.epsilon * this.DXX * (ys * as(:)); % \sum_{j=1}^n a(j) D_xx  y(:,j)
            
            % conservative nonlinearity (u^2/2)_x -> D_x * diag(u_exp) * u_imp
            if(this.conservative == 1)
                c_adv = ys * as(:); % vector \sum_{j=1}^n a(j) y(:,j)
                c_adv = spdiags(c_adv, 0, dim, dim);
                
                y = (this.I - this.DX * c_adv / 2) \ (x + c_diff);
            else % nonconservative nonlinearity u u_x -> diag(u_exp) * D_x * u_imp
                c_adv = ys * as(:); % vector \sum_{j=1}^n a(j) y(:,j)
                c_adv = spdiags(this.DX * c_adv, 0, dim, dim);
                
                y = (this.I - c_adv) \ (x + c_diff);
            end
            
        end
        
        function reset(this)
            this.setInitialCondition();
            this.setOperators();
        end
        
        % =================================================================
        % Start ARG LOCK functions
        % =================================================================
        
        function clearLockedArgs(this, num_args_prealloc)
            if(nargin < 2)
                num_args_prealloc = 1;
            end
            
            this.locked_arg_count = 0;
            this.locked_args      = zeros(this.Nx, 2, num_args_prealloc);
            this.locked_arg_inds  = zeros(num_args_prealloc, 1);
        end
        
        function id = lockArg(this, arg, arg_index)
            
            switch(arg_index)                
                case 1                    
                    id = this.lockArg1(arg);
                case 2                    
                    id = this.lockArg2(arg);
            end

        end

        function id = lockArg1(this, arg)
            
            id = this.locked_arg_count + 1;
            this.locked_arg_inds(id) = 1;
            
            if(this.conservative == 1) % -- conservative F ----------------  
                this.locked_args(:, 1, id) = this.epsilon * this.DXX * arg;
                this.locked_args(:, 2, id) = arg;          
            else % ------------------------ non conservative F ------------  
                this.locked_args(:, 1, id) = this.epsilon * this.DXX * arg;
                this.locked_args(:, 2, id) = this.DX * arg;
            end
            
            this.locked_arg_count = this.locked_arg_count + 1;

        end

        function lockArg2(this, args)
            
            error('not implemented');

        end

        function rhs = F_lockedArg(this, y, id)
            
            switch(this.locked_arg_inds(id))                
                case 1                    
                    rhs = this.F_lockedArg1(y, id);
                case 2                    
                    rhs = this.F_lockedArg2(y, id);
            end

        end

        function rhs = F_lockedArg1(this, y_exp, id)
            
            %F ODE right-hand-side
            dxx_y_imp = this.locked_args(:, 1, id);
                
            if(this.conservative == 1)
                y_imp = squeeze(this.locked_args(:, 2, id));
                rhs = dxx_y_imp + this.DX * (y_exp .* y_imp) / 2;
            else
                dx_y_imp = squeeze(this.locked_args(:, 2, id));
                rhs = dxx_y_imp + y_exp .* dx_y_imp;
            end

        end

        function F_lockedArg2(this, y_imp, id)
            
            error('not implemented');

        end

        % =================================================================
        % End ARG LOCK functions
        % =================================================================
        
    end
    
    methods(Access = protected)
       
        function setInitialCondition(this)
            
            f = @(x) exp(-3*x.^2);
            xs = linspace(this.xspan(1), this.xspan(end), this.Nx+2)'; 
            xs([1,end]) = [];
            
            this.xs = xs;
            this.initial_condition = f(xs);
            
        end
        
        function setOperators(this)
            
            Nx = this.Nx;
            dx = diff(this.xspan) / (Nx + 1);
            
            e = ones(Nx,1);
            this.DX = spdiags([-e, e], [-1, 1], Nx, Nx) / (2 * dx);
            this.DXX = spdiags([e, -2*e, e], [-1, 0, 1], Nx, Nx) / (dx^2);
            this.I = speye(Nx);
            
        end
    
    end
    
end