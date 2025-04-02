classdef MR_NPRK_WRAP < Integrator
    %NPRK_IMEX_WRAP Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % reduced underlying explicit integrator (m steps of user provided tableau)
        Ah2
        bh2
        ch2
        % reduced explicit integrator (single step of user provided tableau)
        A2_st
        b2_st
        c2_st
        m2
    end

    methods
        function this = MR_NPRK_WRAP(methodHandle, m)
            %EULER Construct an instance of this class
            %   Detailed explanation goes here
            
            % store reduced tableau
            [this.A2_st, this.b2_st, this.c2_st] = methodHandle();
            this.m2 = m;
            
            % store reduced explicit method repeated tableau
            [this.Ah2, this.bh2, this.ch2] = nStepTableau(methodHandle, m);

        end
        
    end

    methods(Static)

         function [ store_inds_cell ] = generateStoreIndsCell(A, m, store_inds)
         % determines stages indices to save during each step if the RK
         % method (A,b) is repeated m times

            nstages = size(A, 1);
            store_inds_cell = cell(m, 1);
            
            % prepare store inds
            for i = 1 : length(store_inds)
                step_index = ceil(store_inds(i) / nstages);
                stage_index = 1 + mod(store_inds(i) - 1, nstages);
                store_inds_cell{step_index} = [store_inds_cell{step_index}, stage_index];
            end      
            
         end

         function [y0, Ys] = rkSteps(G, A, b, m, y0, h, store_inds_cell)
        % steps using m steps of the RK method A, b

            Ys = cell(1, m);
            for i = 1 : m
                [y0, Ys{i}] = MR_NPRK_WRAP.rkStep(G, A, b, y0, h / m, store_inds_cell{i});
            end
            Ys = cell2mat(Ys);

        end

        function [y0, Ys] = rkStep(G, A, b, y0, h, store_inds)
            % RK implements a generic explicit Runge Kutta integrator of the form:
            %
            %   Y_1 = G(y_n)
            %   Y_i = y_n + \sum_{j=1}^{i-1} h a_{ij} G(Y_j)    i = 2, ..., s
            %
            %   y_{n+1} = y_n + \sum_{j=1}^s h b_j G(Y_j)
            %
            % PARAMETERS
            %   G - rhs function
            %   A - RK matrix
            %   b - RK weights
            %   w - initial condition
            %   h - stepsize
            %   store_inds - array of integers containing indices of stages
            %                that should be stored

            nstages = size(A, 1);
            K  = zeros(length(y0), nstages);
            Ys = zeros(length(y0), length(store_inds));
            store_flag  = false(nstages, 1); store_flag(store_inds) = true; 
            store_count = 1;

            % -- compute stage values -------------------------------------
            K(:,1) = G(y0);
            for i = 2 : nstages
                K(:, i) = y0;
                for j = 1 : i - 1
                    K(:, i) = K(:, i) + h * A(i, j) * K(:, j);
                end
                if(store_flag(i)) % store specified stage values
                    Ys(:, store_count) = K(:, i);
                    store_count = store_count + 1;
                end
                K(:, i) = G(K(:, i));
            end

            % -- compute output -------------------------------------------
            for i = 1 : nstages
                y0 = y0 + h * b(i) .* K(:, i);
            end

        end

    end

end
