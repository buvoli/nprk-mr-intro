classdef WENOFD5 < handle
    %WENOFD5 fifth-order finite difference weno method (periodic)
    
    properties

        Nx  % number of grid points
        dx  % grid spacing (assumes uniform)
        W5M % weno matrices

    end
    
    methods
        
        function obj = WENOFD5(Nx, dx)
            %WENOTEST Construct an instance of this class
            %   Nx - number of spatial grid points
            %   dx - grid spacing (assumed uniform)
            obj.Nx  = Nx;
            obj.dx  = dx;
            obj.W5M = WENOFD5.WENO5Matrices(Nx);
        end
        
        function u_x_op = DX_upwindOperator(this, u)
            % returns du/dx using upwinding weno
            
            % evaluate smoothness indicator stencils
            SI1 = this.W5M{1,1} * u;
            SI2 = this.W5M{1,2} * u;
            SI3 = this.W5M{1,3} * u;
            SI4 = this.W5M{1,4} * u;
            
            % evaluate derivative stencils
            DS1_op = this.W5M{2,1};
            DS2_op = this.W5M{2,2};
            DS3_op = this.W5M{2,3};

            u_x_op = WENOFD5.computeUpwindDx(SI1, SI2, SI3, SI4, DS1_op, DS2_op, DS3_op, this.Nx, this.dx);
            %u_x_op = WENOFD5.UpwindDxOperator(SI1, SI2, SI3, SI4, DS1_op, DS2_op, DS3_op, this.Nx, this.dx);       
        
        end

        function u_x_op = DX_downwindOperator(this, u)
            % returns du/dx using upwinding weno
            
            % evaluate smoothness indicator stencils
            SI1 = this.W5M{1,1} * u;
            SI2 = this.W5M{1,2} * u;
            SI3 = this.W5M{1,3} * u;
            SI4 = this.W5M{1,4} * u;
            
            % evaluate derivative stencils
            DS2_op = this.W5M{2,2};
            DS3_op = this.W5M{2,3};
            DS4_op = this.W5M{2,4};

            u_x_op = WENOFD5.computeDownwindDx(SI1, SI2, SI3, SI4, DS2_op, DS3_op, DS4_op, this.Nx, this.dx); 
            %u_x_op = WENOFD5.DownwindDxOperator(SI1, SI2, SI3, SI4, DS2_op, DS3_op, DS4_op, this.Nx, this.dx);
            
        end

        function u_x_op = gDxOperator(this, u, g)
            % returns g * du/dx and applys upwinding if g >= 0 and 
            % downwinding if g < 0

            upwind_op = this.DX_upwindOperator(u);
            downwind_op = this.DX_downwindOperator(u);

            upwind_mask = spdiags(g >= 0, 0, this.Nx, this.Nx);
            downwind_mask = spdiags(g < 0, 0, this.Nx, this.Nx);

            u_x_op = upwind_mask * upwind_op + downwind_mask * downwind_op;

        end

        function u_x = DX_upwind(this, u)
            % returns du/dx using upwinding weno
            
            % evaluate smoothness indicator stencils
            SI1 = this.W5M{1,1} * u;
            SI2 = this.W5M{1,2} * u;
            SI3 = this.W5M{1,3} * u;
            SI4 = this.W5M{1,4} * u;
            
            % evaluate derivative stencils
            DS1 = this.W5M{2,1} * u;
            DS2 = this.W5M{2,2} * u;
            DS3 = this.W5M{2,3} * u;

            u_x = WENOFD5.computeUpwindDx(SI1, SI2, SI3, SI4, DS1, DS2, DS3, this.Nx, this.dx);       
        
        end

        function u_x = DX_downwind(this, u)
            % returns du/dx using downwinding weno

            % evaluate smoothness indicator stencils
            SI1 = this.W5M{1,1} * u;
            SI2 = this.W5M{1,2} * u;
            SI3 = this.W5M{1,3} * u;
            SI4 = this.W5M{1,4} * u;
            
            % evaluate derivative stencils
            DS2 = this.W5M{2,2} * u;
            DS3 = this.W5M{2,3} * u;
            DS4 = this.W5M{2,4} * u;

            u_x = WENOFD5.computeDownwindDx(SI1, SI2, SI3, SI4, DS2, DS3, DS4, this.Nx, this.dx);

        end
          
        function u_x = gDx(this, u, g)
            % returns g * du/dx and applys upwinding if g >= 0 and 
            % downwinding if g < 0
            
            % evaluate smoothness indicator stencils
            SI1 = this.W5M{1,1} * u;
            SI2 = this.W5M{1,2} * u;
            SI3 = this.W5M{1,3} * u;
            SI4 = this.W5M{1,4} * u;
            
            % evaluate derivative stencils
            DS1 = this.W5M{2,1} * u;
            DS2 = this.W5M{2,2} * u;
            DS3 = this.W5M{2,3} * u;
            DS4 = this.W5M{2,4} * u;
            
            ux_up = WENOFD5.computeUpwindDx(SI1, SI2, SI3, SI4, DS1, DS2, DS3, this.Nx, this.dx);
            ux_dw = WENOFD5.computeDownwindDx(SI1, SI2, SI3, SI4, DS2, DS3, DS4, this.Nx, this.dx);
            u_x   = (g >= 0) .* ux_up + (g < 0) .* ux_dw;
        
        end

    end

    methods (Static)
            
        function u_x = computeUpwindDx(SI1, SI2, SI3, SI4, DS1, DS2, DS3, Nx, dx) % upwind WENO
        
            B0 = WENOFD5.B( 13 / 12, SI1, -1, (1/4), SI2, -1, Nx);
            B1 = WENOFD5.B( 13 / 12, SI1,  0, (1/4), SI4,  0, Nx);
            B2 = WENOFD5.B( 13 / 12, SI1,  1, (1/4), SI3,  1, Nx);
        
            [w0, w1, w2] = WENOFD5.W(B0, B1, B2);
        
            fp_jp = WENOFD5.F(w0, w1, w2, DS1, -1, DS2, 0, DS3, 1, Nx);
            fp_jm = fp_jp(WENOFD5.pInds(Nx, -1), :);
        
            u_x = (fp_jp - fp_jm) / dx;
        
        end
        
        function u_x = computeDownwindDx(SI1, SI2, SI3, SI4, DS2, DS3, DS4, Nx, dx) % downwind WENO
        
            B0 = WENOFD5.B( 13 / 12, SI1, 2, (1/4), SI3,  2, Nx);
            B1 = WENOFD5.B( 13 / 12, SI1, 1, (1/4), SI4,  1, Nx);
            B2 = WENOFD5.B( 13 / 12, SI1, 0, (1/4), SI2,  0, Nx);
        
            [w0, w1, w2] = WENOFD5.W(B0, B1, B2);
        
            fm_jp = WENOFD5.F(w2, w1, w0, DS2, 0, DS3, 1, DS4, 2, Nx);
            fm_jm = fm_jp(WENOFD5.pInds(Nx, -1), :);
            
            u_x = (fm_jp - fm_jm) / dx;
        
        end

        function b = B(c1, S1, shift1, c2, S2, shift2, N)
        % computes smoothness indicators   
            ind1 = WENOFD5.pInds(N, shift1);
            ind2 = WENOFD5.pInds(N, shift2);
            b = c1 * S1(ind1, :) .^2 + c2 * S2(ind2, :) .^2;
        end

        function [w0,w1,w2] = W(B0, B1, B2) 
        % computes weno weights
        
            % parameters
            gamma0 = 1 / 10;
            gamma1 = 6 / 10;
            gamma2 = 3 / 10;
            epsilon = 1e-6;
        
            A0 = gamma0 ./ (epsilon + B0) .^ 2;
            A1 = gamma1 ./ (epsilon + B1) .^ 2;
            A2 = gamma2 ./ (epsilon + B2) .^ 2;
            Atot = (A0 + A1 + A2);
            
            w0 = A0 ./ Atot;
            w1 = A1 ./ Atot;
            w2 = A2 ./ Atot;
            
        end

        function f = F(W0, W1, W2, S1, shift1, S2, shift2, S3, shift3, N)
        % computes flux 
            ind1 = WENOFD5.pInds(N, shift1);
            ind2 = WENOFD5.pInds(N, shift2);
            ind3 = WENOFD5.pInds(N, shift3);
            f = W0 .* S1(ind1, :) + W1 .* S2(ind2, :) + W2 .* S3(ind3, :);
        end

        function i = pInds(N, shift)
        % periodic indices
            
            % inefficient simple code
            %i = mod((1:N) + shift - 1, N) + 1;
        
            % efficient code (valid for |shift| < N)
            if(shift >= 0)
                i = [shift + 1 : N, 1 : shift];
            else
                i = [N + shift + 1 : N, 1 : N + shift];
            end

        end

        function W5M = WENO5Matrices(Nx)
        % Weno5 derivative operator for d/dx on periodic domain
            
            % Smoothness Indicator Stencils
            SS1 = WENOFD5.sTriCirc([ 1, -2,  1], Nx);
            SS2 = WENOFD5.sTriCirc([ 1, -4,  3], Nx);
            SS3 = WENOFD5.sTriCirc([ 3,  -4,  1], Nx);
            SS4 = WENOFD5.sTriCirc([ 1,  0, -1], Nx);
            
            % Derivative Stencils
            DS1 = WENOFD5.sTriCirc([ 2,  -7,  11] / 6, Nx);
            DS2 = WENOFD5.sTriCirc([-1,   5,  2 ] / 6, Nx);
            DS3 = WENOFD5.sTriCirc([ 2,   5,  -1] / 6, Nx);
            DS4 = WENOFD5.sTriCirc([ 11, -7,  2 ] / 6, Nx);
            
            W5M = {
                SS1, SS2, SS3, SS4;
                DS1, DS2, DS3, DS4
            };
        
        end

        function TC = sTriCirc(w, N, shift)
        % returns a sparse, tridiagonal circulant matrix
        % if shift is specified, it is equivalent to doing a cylical
        % permutation of the spatial variables 
        
            if(nargin < 3)
                shift = 0;
            end

            if(shift == 0)
                TC  = spdiags([w(3), w(1), w(2), w(3), w(1)], [-N+1, -1:1, N-1], N, N);
            else
                TC  = spdiags(repmat([w(1), w(2), w(3)], 1, 3), [-(N+1):-N+1 -1:1, N-1:(N+1)] + shift, N, N);
            end
            
        end

        %  === Start Operator Functions ===================================

        function P = cyclicPermutationMatrix(shift, N)
            if(shift > 0)
                P = spdiags( [ 1, 1 ], [ - N + shift, shift ], N, N );
            elseif(shift < 0)
                P = spdiags( [ 1, 1 ], [ shift, N + shift ], N, N );
            else
                P = speye(N);
            end        
        end
        
        function u_x = UpwindDxOperator(SI1, SI2, SI3, SI4, DS1_op, DS2_op, DS3_op, Nx, dx) % upwind WENO
        
            B0 = WENOFD5.B( 13 / 12, SI1, -1, (1/4), SI2, -1, Nx);
            B1 = WENOFD5.B( 13 / 12, SI1,  0, (1/4), SI4,  0, Nx);
            B2 = WENOFD5.B( 13 / 12, SI1,  1, (1/4), SI3,  1, Nx);
        
            [w0, w1, w2] = WENOFD5.W(B0, B1, B2);
        
            fp_jp = WENOFD5.FOperator(w0, w1, w2, DS1_op, -1, DS2_op, 0, DS3_op, 1, Nx);
            fp_jm = WENOFD5.cyclicPermutationMatrix(-1, Nx) * fp_jp;
        
            u_x = (fp_jp - fp_jm) / dx;
        
        end

        function u_x = DownwindDxOperator(SI1, SI2, SI3, SI4, DS2_op, DS3_op, DS4_op, Nx, dx) % downwind WENO
        
            B0 = WENOFD5.B( 13 / 12, SI1, 2, (1/4), SI3,  2, Nx);
            B1 = WENOFD5.B( 13 / 12, SI1, 1, (1/4), SI4,  1, Nx);
            B2 = WENOFD5.B( 13 / 12, SI1, 0, (1/4), SI2,  0, Nx);
        
            [w0, w1, w2] = WENOFD5.W(B0, B1, B2);
        
            fm_jp = WENOFD5.F(w2, w1, w0, DS2_op, 0, DS3_op, 1, DS4_op, 2, Nx);
            fm_jm = WENOFD5.cyclicPermutationMatrix(-1, Nx) * fm_jp;
            
            u_x = (fm_jp - fm_jm) / dx;
        
        end

        function f_op = FOperator(W0, W1, W2, S1, shift1, S2, shift2, S3, shift3, N)
        % computes flux operator
            
            P1 = WENOFD5.cyclicPermutationMatrix(shift1, N);
            P2 = WENOFD5.cyclicPermutationMatrix(shift2, N);
            P3 = WENOFD5.cyclicPermutationMatrix(shift3, N);
        
            f_op = spdiags(W0, 0, N, N) * S1 * P1 + ...
                   spdiags(W1, 0, N, N) * S2 * P2 + ...
                   spdiags(W2, 0, N, N) * S3 * P3;
            
        end

    end

end
