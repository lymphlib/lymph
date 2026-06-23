%> @file  PicardIterations.m
%> @author Mattia Corti
%> @date 11 May 2026
%> @brief Implementation of Picard iterations for the Fisher-KPP problem.
%>
%==========================================================================
%> @section classPicardIterationsFKPP Class description
%==========================================================================
%> @brief            Implementation of Picard iterations for the resolution
%> of the Fisher-KPP problem (the solver adapts to semi-implicit choice).
%>
%> @param Data              Struct with problem's data
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param Matrices          Matrices associated to the linear terms
%> @param Solution          Solution at the previous timesteps
%> @param F_old             Forcing term associated to the previous timestep
%> @param F_new             Forcing term associated to the current timestep
%>
%> @retval c_star	    Solution of the Picard iterations
%>
%==========================================================================

function [c_star] = PicardIterations(Data, femregion, Matrices, Solution, F_old, F_new)

    % Initialization of the Picard iterations
    ii  = 0;
    err = Data.NLS_tolerance + 1;

    switch Data.NonLinearSolver 
        case 'Semi-implicit'
            c_star = Solution.c_old + 0.5*(Data.theta==0.5)*(Solution.c_old-Solution.c_oold);
            NL = 0;
        case 'Picard iterations'
            c_star = Solution.c_old;
            NL = 1;
            % Nonlinear matrix term associated to the previous timestep
    end
    
    [MN_c_old] = ContractTrilinearMatrix(Matrices.NonLinear.M_NL, c_star, femregion.nel, femregion.nbases);
        
    % Picard iterations scheme
    while (NL || ~ii) && (err > Data.NLS_tolerance) && (ii < Data.NLS_max_it)
        
        % Nonlinear matrix term associated to the current timestep
        [MN_c_star] = ContractTrilinearMatrix(Matrices.NonLinear.M_NL, c_star, femregion.nel, femregion.nbases); 
    
        % Construction of LHS and RHS
        LHS =  Matrices.M_prj + Data.dt*Data.theta*Matrices.A - Data.theta*Data.dt*(Matrices.M-(1+NL*(Data.theta-1))*MN_c_star-(Data.theta==0.5)*NL*MN_c_old);

        F   = (Matrices.M_prj - Data.dt*(1-Data.theta)*Matrices.A + (1-Data.theta)*Data.dt*(Matrices.M-(1-0.5*NL)*MN_c_old))*Solution.c_old ...
            + Data.dt*Data.theta*F_new.F + Data.dt*(1-Data.theta)*F_old.F;
    
        % Resolution of the linear system
        c_h = LHS\F;
    
        % Fixed point iteration update
        ii = ii + 1;
        
        % Solution relaxation
        c_app = (1-NL)*c_h+NL*(Data.NLS_relax*c_h + NL*(1-Data.NLS_relax)*c_star);

        % Error computation
        err = (c_star-c_app)'*Matrices.dGA*(c_star-c_app)/(c_star'*Matrices.dGA*c_star);
        
        % Solution update
        c_star = c_app;

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
