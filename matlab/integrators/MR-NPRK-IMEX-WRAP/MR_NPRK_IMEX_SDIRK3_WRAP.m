classdef MR_NPRK_IMEX_SDIRK3_WRAP < MR_NPRK_WRAP
    %NPRK_IMEX_SDIRK3_WRAP third order IMEX MR-NPRK implicity wrapped 
    %method using three stage L-stable S-dirk method

    properties
        name = 'IMEX-NPRK-SDIRK3-WRAP'
    end

    properties(SetAccess = protected)
        % method configuration
        omega
        variant
        % method coefficients
        coupling_coeff
        % reduced underlying implicit integrator
        Ah1
        bh1
        ch1
        % misc
        store_inds_cell % stores stages to save during explicit timestepping
    end

    methods
        
        function this = MR_NPRK_IMEX_SDIRK3_WRAP(methodHandle, m, variant, omega, name)
            %NPRK_IMEX_SDIRK2_WRAP Construct an instance of this class
            %   Detailed explanation goes here
            
            if(nargin < 3)
                variant = 1;
            end
            if(nargin < 4)
                omega = 2;
            end

            this = this@MR_NPRK_WRAP(methodHandle, m);
            
            [this.Ah1, this.bh1, this.ch1] = this.SDIRK3_Sa();
            this.variant = variant;
            this.omega = omega;
            
            s2 = size(this.A2_st, 1) * this.m2; % num explicit stages
            
            if(this.variant == 1)
                this.coupling_coeff = MR_NPRK_IMEX_SDIRK3_WRAP.couplingCoeffV1(this.bh2, this.ch2, this.omega);
                store_inds = [ s2 - 1, s2 ];
            elseif(this.variant == 2)
                this.coupling_coeff = MR_NPRK_IMEX_SDIRK3_WRAP.couplingCoeffV2(this.bh2, this.ch2, this.omega);
                store_inds = [ s2 ];
            else
                error('invalid variant');
            end

            % ---> set store inds
            this.store_inds_cell = this.generateStoreIndsCell(this.A2_st, this.m2, store_inds);

            if(nargin == 5) % optional name override
                this.name = name;
            end
            
        end

        function ynp = step(this, yn, h, problem)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            switch(this.variant)
                case 1
                    ynp = stepV1(this, yn, h, problem);
                case 2
                    ynp = stepV2(this, yn, h, problem);
            end

        end

        function ynp = stepV1(this, yn, h, problem)

            % ----> local variable copies
            Ah1 = this.Ah1; 
            bh2 = this.bh2;

            a_322         = this.coupling_coeff.a_322;
            a_s2p2_3_s2p1 = this.coupling_coeff.a_s2p2_3_s2p1;
            a_s_3_1       = this.coupling_coeff.a_s_3_1;
            a_s_3_sm1     = this.coupling_coeff.a_s_3_sm1;

            % ----> first implicit stage
            b  = yn;
            Y2 = problem.invertFComp(b, h * Ah1(1,1), yn, yn, 1);

            % ----> second implicit stage
            F_Y2_yn = problem.F(Y2,yn);
            F_Y2_Y2 = problem.F(Y2,Y2);
            b  = yn + h * ( (Ah1(2,1) - a_322) * F_Y2_yn + a_322 * F_Y2_Y2 );
            if(this.omega == 1)
                ys = yn;
            else
                ys = Y2;
            end            
            Y3 = problem.invertFComp(b, h * Ah1(2,2), ys, yn, 1);

            % ----> explicit stages
            if(problem.arg_lock_enabled)
                problem.clearLockedArgs();
                problem.lockArg(Y2, 1);
                G = @(w) problem.F_lockedArg(w, 1);
            else
                G = @(w) problem.F(Y2, w);
            end
            [w1, W] = this.rkSteps(G, this.A2_st, this.b2_st, this.m2, yn, h, this.store_inds_cell);

            % ----> third implicit stage

            F_Y2_Ws2m1 = problem.F(Y2, W(:,1)); % could be recycled from RKsteps
            F_Y3_Ws2m1 = problem.F(Y3, W(:,1));

            Y_sm1 = W(:,2) + h * a_s2p2_3_s2p1 * (F_Y3_Ws2m1 - F_Y2_Ws2m1); 
            
            F_Y2_Ws2  = problem.F(Y2, W(:,2));  % could be recycled from RKsteps
            F_Y2_Ysm1 = problem.F(Y2, Y_sm1);
            F_Y3_yn   = problem.F(Y3, yn);
            F_Y3_Ysm1 = problem.F(Y3, Y_sm1);

            b = w1 ...
                - h * bh2(end) * F_Y2_Ws2 ...
                + h * a_s_3_1 * (F_Y3_yn - F_Y2_yn) ...
                + h * a_s_3_sm1 * F_Y3_Ysm1 ...
                + h * (bh2(end) - a_s_3_sm1 - Ah1(3,3)) * F_Y2_Ysm1;
            
            % -------------------------------------------------------------
            %  Easier to read definition of b (slightly less efficient)
            % -------------------------------------------------------------
            % b = w1 ...
            %     + h * bh2(end) * (F_Y2_Ysm1 - F_Y2_Ws2) ...
            %     + h * a_s_3_1 * (F_Y3_yn - F_Y2_yn) ...
            %     + h * a_s_3_sm1 * (F_Y3_Ysm1 - F_Y2_Ysm1) ...
            %     - h * Ah1(3,3) * F_Y2_Ysm1;
            % -------------------------------------------------------------

            ynp = problem.invertFComp(b, h * Ah1(3,3), Y_sm1, yn, 1);

            % -------------------------------------------------------------
            %  More readable code for third implicit stage 
            %  (but its wastes function evaluations)
            % -------------------------------------------------------------
            % Y_sm1 = W(:,2) + h * a_s2p2_3_s2p1 * (problem.F(Y3, W(:,1)) - problem.F(Y2, W(:,1))); 
            % b     = w1 ...
            %     + h * bh2(end) * (problem.F(Y2, Y_sm1) - problem.F(Y2, W(:,2))) ...
            %     + h * a_s_3_1 * (problem.F(Y3, yn) - problem.F(Y2, yn)) ...
            %     + h * a_s_3_sm1 * (problem.F(Y3, Y_sm1) - problem.F(Y2, Y_sm1)) ...
            %     - h * Ah1(3,3) * problem.F(Y2, Y_sm1);
            % 
            % ynp = problem.invertFComp(b, h * Ah1(3,3), Y_sm1, yn, 1);
            % -------------------------------------------------------------

        end

        function ynp = stepV2(this, yn, h, problem)

            % ----> local variable copies
            Ah1 = this.Ah1; 
            bh2 = this.bh2;
            
            a_322     = this.coupling_coeff.a_322;
            delta     = this.coupling_coeff.delta;
            a_s_3_1   = this.coupling_coeff.a_s_3_1;
            a_s_3_sm1 = this.coupling_coeff.a_s_3_sm1;
            
            % ----> first implicit stage
            b  = yn;
            Y2 = problem.invertFComp(b, h * Ah1(1,1), yn, yn, 1);
            
            % ----> second implicit stage
            b  = yn + h * ( (Ah1(2,1) - a_322) * problem.F(Y2,yn) + a_322 * problem.F(Y2,Y2) );
            if(this.omega == 1)
                ys = yn;
            else
                ys = Y2;
            end            
            Y3 = problem.invertFComp(b, h * Ah1(2,2), ys, yn, 1);

            % ----> explicit stages
            if(problem.arg_lock_enabled)
                problem.clearLockedArgs();
                id_2 = problem.lockArg(Y2, 1);
                id_3 = problem.lockArg(Y3, 1);
                G = @(w) (1 - delta) * problem.F_lockedArg(w, id_2) + delta * problem.F_lockedArg(w, id_3);
            else
                G = @(w) (1 - delta) * problem.F(Y2, w) + delta * problem.F(Y3, w);
            end
            [w1, W] = this.rkSteps(G, this.A2_st, this.b2_st, this.m2, yn, h, this.store_inds_cell);
            
            % ----> third implicit stage
            F_Y2_W = problem.F(Y2, W);
            H_yn   = (problem.F(Y3, yn) - problem.F(Y2, yn));
            H_s2   = (problem.F(Y3, W) - F_Y2_W);
            
            cr  = (a_s_3_1 - delta * bh2(1)) * H_yn ...
                  + (a_s_3_sm1 - delta * bh2(end)) * H_s2 ...
                  - Ah1(3,3)  * F_Y2_W;           
            b   = w1 + h * cr; 
            ynp = problem.invertFComp(b, h * Ah1(3,3), W, yn, 1);

        end

    end

    methods(Static, Access = protected)

        function [A, b, c] = SDIRK3_Sa()
        %SDIRK3_Sa Stiffly accurate, L-stable three stage SDIRK method
        
        lm = 0.43586652150845899941601945119356;
        
        A = [
            lm                                  0                               0;
            (1 - lm) / 2                        lm                              0;
            (-6 * lm^2 + 16 * lm - 1) / 4       (6 * lm^2 - 20 * lm + 5) / 4    lm;
            ];
        
        b = A(3,:);
        
        c = sum(A, 2);
        
        end

        function [coeff] = couplingCoeffV1(bh2, ch2, omega)
        % coupling coefficients for MR-NPRK variant one (Table X) in paper
        %
        % Parameters:
        % bh2 - b vector for reduced underlying explicit integrator
        % ch2 - c vector for reduced underlying explicit integrator
        % omega - second implicit solve is Y_3 = b + F(Y_3, Y_omega)
        
        [Ah1, ~, ch1] = MR_NPRK_IMEX_SDIRK3_WRAP.SDIRK3_Sa();

        % ----> coefficient a_{3,2,2}
        if(omega == 1)
            a_322 = 0.1825367204468751798095174531383;
        elseif(omega == 2)
            a_322 = - 0.253329801061583819606501998055;
        else
            error('invalid omega');
        end

        % ----> coefficient a_{\scmp{2}+2,3,\scmp{2}+1}
        a_s2p2_3_s2p1 = (1 - 3 * ch1(1)) / ( 6 * bh2(end) * (ch1(2) - ch1(1)));
        
        % ----> coefficient a_{s,3,1}
        n = 1 / 3 - Ah1(3,3) * ch2(end) + ch1(1) * (Ah1(3,3) * ch2(end) - 1 / 2) - Ah1(3,2) * (ch1(2) - ch1(1)) * ch2(end);
        d = ch2(end) * (ch1(1) - ch1(2));
        a_s_3_1 = n / d;
        
        % ----> coefficient a_{s,3,s-1}
        a_s_3_sm1 = Ah1(3,2) - a_s_3_1;

        % ----> coefficient struct
        coeff = struct(                         ...
            'a_322',            a_322,           ...
            'a_s2p2_3_s2p1',    a_s2p2_3_s2p1,  ...
            'a_s_3_1',          a_s_3_1,        ...
            'a_s_3_sm1',        a_s_3_sm1       ...
        );

        end

        function [ coeff ] = couplingCoeffV2(bh2, ch2, omega)
        % coupling coefficients for MR-NPRK variant two (Table X) in paper
        %
        % Parameters:
        % bh2 - b vector for reduced underlying explicit integrator
        % ch2 - c vector for reduced underlying explicit integrator
        % omega - second implicit solve is Y_3 = b + F(Y_3, Y_omega)
        
        [Ah1, ~, ch1] = MR_NPRK_IMEX_SDIRK3_WRAP.SDIRK3_Sa();

        % ----> coefficient a_{3,2,2}
        if(omega == 1)
            a_322 = 0.1825367204468751798095174531383;
        elseif(omega == 2)
            a_322 = - 0.253329801061583819606501998055;
        else
            error('invalid omega');
        end

        % ----> coefficient delta
        delta = - 0.36350683689006809456061097715969; % = (1 - 3 * c1(1))/(3 * (c1(2) - c1(1)));
        
        % ----> coefficient a_{s,3,1}
        X1 = Ah1(3,2) - delta * (1 - bh2(1) - bh2(end));
        X3 = bh2(end) - X1;
        X4 = (1/2) - (bh2(end) * ch2(end));

        N1 = ((1 - delta) * X4 + ch2(end) * X3 - Ah1(3,3) * ch2(end));
        N2 = (delta * X4 + ch2(end) * X1);

        num = 1 / 3 - ch2(end) * Ah1(3,3) - ch1(1) * N1 - ch1(2) * N2;
        den = (ch1(1) - ch1(2)) * (ch2(end));

        a_s_3_1 = num / den;
        
        % ----> coefficient a_{s,3,s-1}
        a_s_3_sm1 = Ah1(3, 2) - a_s_3_1 - delta * (1 - bh2(1) - bh2(end));

        % ----> coefficient struct
        coeff = struct(                 ...
            'a_322',        a_322,      ...
            'delta',        delta,      ...
            'a_s_3_1',      a_s_3_1,    ...
            'a_s_3_sm1',    a_s_3_sm1   ...
        );

        end

    end

end
