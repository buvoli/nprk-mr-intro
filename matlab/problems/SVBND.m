classdef SVBND < Problem
    %SVAD spatially varying advection-diffusion equation
    %   u_t + nu * a(x) u_x + epsilon * b(x) u_xx
    
    properties
        % Generic properties
        tspan = [0 5]
        initial_condition
        % Spatial grid
        Nx    = 300         % number of internal gridpoints
        xspan = [-2, 2]     % spatial domain range
        xs                  % spatial grid
        % PDE parameters
        gamma = 1 / 2;
        nu = 1;
        epsilon = 1;
        conservative = 1;
        arg_lock_enabled = false;
        initial_condition_type = 1 % 1 - gaussian, 2 - three gaussians, 3 - three gaussians
    end
    
    properties (SetAccess = protected)
        % operators
        NLD
        W5
        epsBx
        nuAx
        I   
    end
    
    methods
        
        function this = SVBND()
            %BURGERS Construct an instance of this class
            %   Detailed explanation goes here
            this.reset();
        end

        function rhs = F(this, y_imp, y_exp)
            %F ODE right-hand-side
            
            rhs_diff = this.epsBx * this.NLD.applyOpPeriodic(abs(y_exp), y_imp);
            if(this.conservative)
                rhs_adv = this.nuAx * this.W5.DX_downwind(y_exp.^2);
            else                
                rhs_adv = this.nuAx * y_imp .* this.W5.DX_downwind(y_exp);
            end            
            rhs = rhs_adv + rhs_diff;
            
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
            %INVERTIMP1 solves the equation y = x + \sum_{j=1}^n a_j H(y,y_j)            
            
            num_as = length(as);
            
            % --> implicit operator
            op_imp = this.I;
            for i = 1 : num_as
                op_imp = op_imp - as(i) * this.epsBx * this.NLD.OpPeriodic(abs(ys(:,i)));
            end

            % --> explicit operator            
            b_exp = 0;
            for i = 1 : num_as
                b_exp = b_exp + as(i) * this.nuAx * this.W5.DX_downwind(ys(:,i).^2);
            end

            y = op_imp \ (x + b_exp);
            
        end
        
        function y = invertImp2(this, x, as, ys, y0)
            %INVERTIMP2 solves the equation y = x + \sum_{j=1}^n a_j H(y_j,y)            
            
            error('not implemented');
            
        end
        
        function reset(this)
            this.setInitialCondition();
            this.setOperators();
        end
        
    end
    
    methods(Access = protected)
       
        function setInitialCondition(this)
            
            switch(this.initial_condition_type)
                case 1
                    f = @(x) exp(-2 * x.^2);
                    this.tspan = [0 8];
                case 2
                    f = @(x) (1e-2) + exp(-60 * (x + 3/2).^2) + exp(-60 * x.^2);% + exp(-60 * (x - 3/2).^2);
                    this.tspan = [0 5];
                case 3
                    f = @(x) (1e-2) + exp(-60 * (x + 3/2).^2) + exp(-60 * x.^2) + exp(-60 * (x - 3/2).^2);
                    this.tspan = [0 5];
            end

            xs = linspace(this.xspan(1), this.xspan(end), this.Nx+1)'; 
            xs(end) = [];
            
            this.xs = xs;
            this.initial_condition = f(xs);
            
        end
        
        function setOperators(this)
            
            Nx = this.Nx;
            dx = diff(this.xspan) / (Nx + 1);

            this.W5    = WENOFD5(Nx, dx);
            this.nuAx  = spdiags(this.nu * this.a(this.xs), 0, Nx, Nx);
            
            this.NLD   = NLD2(@(u) u.^this.gamma, Nx, dx);
            this.epsBx = spdiags(this.epsilon * this.b(this.xs), 0, Nx, Nx);

            this.I     = speye(Nx);
 
        end

        function ax = a(this, x)
            
            x_center = mean(this.xspan);
            x_width  = diff(this.xspan);

            g_mean = x_center - x_width / 4;
            g_std  = x_width / 30;            
            ax     = 0.5 + 2 * exp( -1 * (x - g_mean).^2 / (2 * g_std^2));

        end

        function bx = b(this, x)
            
            x_center = mean(this.xspan);
            x_width  = diff(this.xspan);

            g_mean = x_center + x_width / 4;
            g_std  = x_width / 30;            
            bx     = exp( -1 * (x - g_mean).^2 / (2 * g_std^2));

            % clip very small diffusion values
            bx(abs(bx) < 1e-10) = 0;

        end
    
    end
    
end