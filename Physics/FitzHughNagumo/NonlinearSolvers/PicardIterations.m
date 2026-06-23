%> @file  PicardIterations.m
%> @author Mattia Corti, Caterina Leimer
%> @date 11 May 2026
%> @brief Implementation of Picard iterations for the FHN problem.
%>
%==========================================================================
%> @section classPicardIterationsFHN Class description
%==========================================================================
%> @brief            Implementation of Picard iterations for the resolution
%> of the FHN problem (the solver adapts to semi-implicit choice).
%>
%> @param Data              Struct with problem's data
%> @param mesh              Mesh struct
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param Matrices          Matrices associated to the linear terms
%> @param Solution_old      Solution at the previous timesteps
%> @param F_old             Forcing term associated to the previous timestep
%> @param F_new             Forcing term associated to the current timestep
%> @param t                 Current time
%>
%> @retval c_star	    Solution of the Picard iterations
%>
%==========================================================================

function [Solution_it,F_new] = PicardIterations(Data, mesh, femregion, Matrices, Solution_old, F_old, F_new, t)

    % Initialization of the Picard iterations
    ii  = 0;
    err = Data.NLS_tolerance + 1;

    Solution_it = Solution_old;
    
    % Picard iterations scheme
    while (strcmp(Data.NonLinearSolver,'Picard iterations') || ~ii) && (err > Data.NLS_tolerance) && (ii < Data.NLS_max_it)
        
        % Construction of LHS and RHS
        LHS_u =  Matrices.M_u + Data.dt*Data.theta*Matrices.A;

        F_u   = (Matrices.M_u - Data.dt*(1-Data.theta)*Matrices.A)*Solution_old.u_h +...
                Data.dt*Data.theta*F_new.Iext + Data.dt*(1-Data.theta)*F_old.Iext + ...
                Data.dt*F_new.F;

        LHS_w = Matrices.M_w;

        F_w = Matrices.M_w*Solution_old.w_h + Data.dt*F_new.G;
    
        % Resolution of the linear system
        Solution.u_h = LHS_u\F_u;
        Solution.w_h = LHS_w\F_w;
    
        % Fixed point iteration update
        ii = ii + 1;
        
        % Solution relaxation
        if strcmp(Data.NonLinearSolver,'Picard iterations')
            
            u_app = (Data.NLS_relax*Solution.u_h + (1-Data.NLS_relax)*Solution_it.u_h);
            w_app = (Data.NLS_relax*Solution.w_h + (1-Data.NLS_relax)*Solution_it.w_h);

            % Error computation
            err = (Solution_it.u_h-u_app)'*Matrices.dGA*(Solution_it.u_h-u_app)/(Solution_it.u_h'*Matrices.dGA*Solution_it.u_h);

            Solution_it.u_h = u_app;
            Solution_it.w_h = w_app;
            
            if Data.theta == 0.5
                Solution_star.u_h = 0.5*Solution_it.u_h+0.5*Solution_old.u_h;
                Solution_star.w_h = 0.5*Solution_it.w_h+0.5*Solution_old.w_h;
            else
                Solution_star = Solution_it;
            end

            [F_new] = ForcingTermAssemblyFHN(Data, mesh, femregion, t, Solution_star, []);

        else

            Solution_it = Solution;

        end

    end

    
    if strcmp(Data.NonLinearSolver,'Picard iterations')
        
        % Display the Picard iterations result
        if err > Data.NLS_tolerance
            fprintf(' - Picard iterations terminated without convergence (final error: %4.3e)\n', err);
        else
            fprintf(' - Picard iterations convergenced in %d iterations (final error: %4.3e)\n', ii, err);
        end
    end

end
