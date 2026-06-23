%> @file  AdaptivityCycle.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 18 February 2026
%> @brief Implementation of adaptivity cycle for the heat equation.
%>
%==========================================================================
%> @section classAdaptivityCycleHeat Class description
%==========================================================================
%> @brief            Implementation of adaptivity cycle for the heat equation.
%>
%> @param Data              Struct with problem's data
%> @param mesh              Mesh struct (neighbor and region structs)
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param Solution	    Structure containing solution and indicators
%> @param Matrices          Matrices associated to the PDE
%> @param F_old             Forcing term associated with the previous timestep
%> @param F_new		    Forcing term associated with the current timestep
%> @param t		    Current time
%>
%> @retval Data		    Struct with problem's data
%> @retval Solution         Struct containing the problem solution
%> @retval femregion        Finite Element struct (see CreateDOF.m)
%> @retval Matrices         Matrices associated to the PDE
%> @retval RHS		    Right-hand-side matrices of the theta-method
%> @retval LHS              Left-hand-side matrices of the theta-method
%> @retval F_new            Forcing term associated with the current timestep
%>
%==========================================================================

function [Data, Solution, femregion, Matrices, RHS, LHS, F_new] = AdaptivityCycle(Data, mesh, femregion, Solution, Matrices, F_old, F_new, t)

    AdaptCount = 0;
    
    while AdaptCount < Data.AdaptIts
    
        fprintf('\nComputing adaptive indicator ... ');
        tic
        %% Indicator computation
        [Solution.Indicator] = ComputeAdaptiveIndicatorHeat(Data, mesh, femregion, Solution, t);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    
        %% Update the polynomial degrees
        [Data, femregion] = ComputeUpdatedDegrees(Data, mesh, femregion, Solution, not(round(t/Data.dt)==Data.AdaptivityStep));
           
        %% Project the solutions on the new femregion
        [Solution] = ProjectSolutionAssemblyHeat(Data, mesh, femregion, t, Solution);
    
        %% Reassemble the Matrices
        fprintf('\nMatrices computation ... \n');
        tic
        [Matrices] = MatrixAssemblyHeat(Data, mesh, femregion, Matrices);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    
        LHS = Matrices.Mprj + Data.dt*Data.theta*(Matrices.A+Matrices.M);
        RHS = Matrices.Mprj - Data.dt*(1-Data.theta)*(Matrices.A+Matrices.M);
    
        %% Right-hand side assembly
        fprintf('\nComputing forcing terms ... \n');
        [F_old] = ForcingTermAssemblyHeat(Data, mesh, femregion, t-Data.dt, F_old);
    
        [F_new] = ForcingTermAssemblyHeat(Data, mesh, femregion, t, F_new);
        
        %% Assembling the dynamic component of the matrices in the two treatments of nonlinear term
    
        % Construction of the complete RHS for the theta method
        F = RHS * Solution.u_old + Data.dt*Data.theta*F_new.F + Data.dt*(1-Data.theta)*F_old.F;
    
        % Linear system resolution
        Solution.u_h = LHS\F;
    
        AdaptCount = AdaptCount + 1;
    end
end
