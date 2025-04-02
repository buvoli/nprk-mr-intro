classdef NLD2
    %NLD2 second-order nonlinear diffusion operator (flux form)
    % efficient and stable descritization of:
    % 
    %                       (k(u) * u_x)_x
    % 
    % Discritization is of the form
    % 
    % (k_{i+1/2}*(u_{i+1} - u_{i}) - k_{i-1/2}*(u_{i} - u_{i-1})) / (h^2)
    % for k_{i+1/2} = (k(u_i) + k(u_{i+1}) / 2
    %
    % This discritzation is equivalent to
    %           
    %   (k(u) * u_x)_x = A(u) * u         A(u) \in \mathbb{R}^{N \times N}
    %
    % -- Dirichlet Boundary Conditions ------------------------------------
    %
    % OpDirichlet(this, u)   -> returns matrix A(u)
    % applyOpDirichlet(u, b) -> matrix free implementation of (k(u) b_x)_x
    %
    % REMARK: use properties drhlt_u_left, drhlt_u_right to set boundaries;
    % default values are zero.
    %
    % -- Periodic Boundary Conditions -------------------------------------
    %
    % OpPeriodic(this, u)   -> returns matrix A(u)
    % applyOpPeriodic(u, b) -> matrix free implementation of (k(u) b_x)_x

    properties

        k   % function handle
        Nx  % number of grid points
        dx  % grid spacing (assumes uniform)

        % boundary conditions (Dirichlet)
        drhlt_u_left  = 0 % left boundary condition
        drhlt_u_right = 0 % right boundary condition

    end
    
    methods
        
        function obj = NLD2(k, Nx, dx)
            %WENOTEST Construct an instance of this class
            %   Nx - number of spatial grid points
            %   dx - grid spacing (assumed uniform)
            
            obj.Nx = Nx;
            obj.dx = dx;
            obj.k  = k;

        end
    
        function Op_u = applyOpDirichlet(this, u, b)
            % computes (k(u) b_x)_x

            tdx2 = (2 * this.dx^2);
            
            ku      = this.k(u) / tdx2;
            ku_left  = this.k(this.drhlt_u_left) / tdx2;
            ku_right = this.k(this.drhlt_u_right) / tdx2;
            ku_mid   = ( ku(1 : end - 1) + ku(2 : end) );
            
            lb_mid = (ku_left  + 2 * ku(1) + ku(2));        % left boundary
            rb_mid = (ku(end-1) + 2 * ku(end) + ku_right);  % right boundary
            diag   = -1 * [lb_mid; ku_mid(1:end-1) + ku_mid(2:end); rb_mid];
            
            Op_u = diag .* b ...                      % diagonal
                   + [0; ku_mid .* b(1 : end - 1)] ... % lower diagonal
                   + [ku_mid .* b(2 : end); 0];        % upper diagonal
        
        end

        function [lower_diag, diag, upper_diag] = OpDirichletDiagonals(this, u)
            % computes diagonsl of operator A such that A * b is (k(u) * b_x)_x
    
            tdx2 = (2 * this.dx^2);
            
            ku       = this.k(u) / tdx2;
            ku_left  = this.k(this.drhlt_u_left) / tdx2;
            ku_right = this.k(this.drhlt_u_right) / tdx2;
            ku_mid   = ( ku(1 : end - 1) + ku(2 : end) );
            
            lb_mid = (ku_left  + 2 * ku(1) + ku(2));        % left boundary
            rb_mid = (ku(end-1) + 2 * ku(end) + ku_right);  % right boundary
            
            lower_diag = [ ku_mid; 0 ];
            diag       = -1 * [lb_mid; ku_mid(1:end-1) + ku_mid(2:end); rb_mid];
            upper_diag = [ 0; ku_mid ];
                    
        end
        
        function Op = OpDirichlet(this, u)
            % computes operator A such that A * b is (k(u) * b_x)_x
    
            tdx2 = (2 * this.dx^2);
            
            ku       = this.k(u) / tdx2;
            ku_left  = this.k(this.drhlt_u_left) / tdx2;
            ku_right = this.k(this.drhlt_u_right) / tdx2;
            ku_mid   = ( ku(1 : end - 1) + ku(2 : end) );
            
            lb_mid = (ku_left  + 2 * ku(1) + ku(2));        % left boundary
            rb_mid = (ku(end-1) + 2 * ku(end) + ku_right);  % right boundary
            
            lower_diag = [ ku_mid; 0 ];
            diag       = -1 * [lb_mid; ku_mid(1:end-1) + ku_mid(2:end); rb_mid];
            upper_diag = [ 0; ku_mid ];
            
            Op = spdiags( ...
                [lower_diag, diag, upper_diag], ...
                [-1, 0, 1], ...
                this.Nx, this.Nx);
        
        end

        function Op_u = applyOpPeriodic(this, u, b)
            % computes (k(u) b_x)_x

            ku       = this.k(u) / (2 * this.dx^2);
            ku_left  = ku(end);
            ku_right = ku(1);
            ku_mid   = ( ku(1 : end - 1) + ku(2 : end) );
            
            lb_mid = (ku_left  + 2 * ku(1) + ku(2));        % left boundary
            rb_mid = (ku(end-1) + 2 * ku(end) + ku_right);  % right boundary
            diag   = -1 * [lb_mid; ku_mid(1:end-1) + ku_mid(2:end); rb_mid];
            
            % periodic corrections
            Op1N = (ku_left + ku(1));
            OpN1 = (ku(end) + ku_right);

            Op_u = diag .* b ...                                    % diagonal
                   + [Op1N * b(end); ku_mid .* b(1 : end - 1)] ...  % lower diagonal
                   + [ku_mid .* b(2 : end); OpN1 * b(1)];           % upper diagonal
        
        end

        function Op = OpPeriodic(this, u)
            % computes operator A such that A * b is (k(u) * b_x)_x
    
            ku       = this.k(u) / (2 * this.dx^2);
            ku_left  = ku(end);
            ku_right = ku(1);
            ku_mid   = ( ku(1 : end - 1) + ku(2 : end) );
            
            lb_mid = (ku_left  + 2 * ku(1) + ku(2));        % left boundary
            rb_mid = (ku(end-1) + 2 * ku(end) + ku_right);  % right boundary
            
            lower_diag = [ ku_mid; 0 ];
            diag       = -1 * [lb_mid; ku_mid(1:end-1) + ku_mid(2:end); rb_mid];
            upper_diag = [ 0; ku_mid ];
            
            Op = spdiags( ...
                [lower_diag, diag, upper_diag], ...
                [-1, 0, 1], ...
                this.Nx, this.Nx);

            Op(1, this.Nx) = (ku_left + ku(1));
            Op(this.Nx, 1) = (ku(end) + ku_right);
        
        end

    end

end

