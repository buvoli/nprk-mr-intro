classdef MR_NPRK_IMEX_SDIRK2_WRAP < MR_NPRK_WRAP
    %NPRK_IMEX_SDIRK2_WRAP second order IMEX MR-NPRK implicity wrapped 
    %method using two stage S-dirk method

    properties
        name = 'IMEX-NPRK-SDIRK2-WRAP'
    end

    properties(SetAccess = protected)
        gamma
        store_inds_cell
    end

    methods
        
        function this = MR_NPRK_IMEX_SDIRK2_WRAP(methodHandle, m, gamma_flag, name)
            %NPRK_IMEX_SDIRK2_WRAP Construct an instance of this class
            %   Detailed explanation goes here
            
            if(nargin < 3)
                gamma_flag = -1;
            end

            this = this@MR_NPRK_WRAP(methodHandle, m);

            % ---> set store inds
            store_inds = [ size(this.A2_st, 1) * this.m2 ];
            this.store_inds_cell = this.generateStoreIndsCell(this.A2_st, this.m2, store_inds);
            
            % ---> set gamma
            if(gamma_flag == 1)
                this.gamma = (2 + sqrt(2)) / 2; % improved stability / damping
            else
                this.gamma = (2 - sqrt(2)) / 2; % improved error
            end

            if(nargin == 4) % optional name override
                this.name = name;
            end
            
        end

        function ynp = step(this, yn, h, problem)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            gm = this.gamma;            
            Y2 = problem.invertFComp(yn, h * gm, yn, yn, 1);
            if(problem.arg_lock_enabled)
                problem.clearLockedArgs();
                problem.lockArg(Y2,1);
                G = @(w) problem.F_lockedArg(w, 1);
            else
                G  = @(w) problem.F(Y2, w);
            end
            [w1, W] = this.rkSteps(G, this.A2_st, this.b2_st, this.m2, yn, h, this.store_inds_cell);
            b = w1 - h * gm * G(W);
            ynp = problem.invertFComp(b, h * gm, W, yn, 1);

        end
    end

end
