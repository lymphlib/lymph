%> @file  AdaptivityCycleFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 18 February 2026
%> @brief Implementation of adaptivity cycle for the FHN equation.
%>
%==========================================================================
%> @section classAdaptivityCycleFHN Class description
%==========================================================================
%> @brief            Implementation of adaptivity cycle for the FHN equation.
%>
%> @param Data              Struct with problem's data
%> @param mesh              Mesh struct (neighbor and region structs)
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param Solution	        Structure containing solution and indicators
%> @param Matrices          Matrices associated to the PDE
%> @param F_old             Forcing term associated with the previous timestep
%> @param F_new		        Forcing term associated with the current timestep
%> @param t		            Current time
%>
%> @retval Data		        Struct with problem's data
%> @retval Solution         Struct containing the problem solution
%> @retval femregion        Finite Element struct (see CreateDOF.m)
%> @retval Matrices         Matrices associated to the PDE
%> @retval F_new            Forcing term associated with the current timestep
%>
%==========================================================================

function [Data, Solution, femregion, Matrices, F_new] = AdaptivityCycleFHN(Data, mesh, femregion, Solution, Matrices, F_old, F_new, t)

    AdaptCount = 0;
    
    while AdaptCount < Data.AdaptIts
    
        fprintf('\nComputing adaptive indicator ... ');
        tic
        %% Indicator computation
        [Solution.t.Indicator] = ComputeAdaptiveIndicatorFHN(Data, mesh, femregion, Solution, t);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    
        %% Update the polynomial degrees
        [Data, femregion] = ComputeUpdatedDegrees(Data, mesh, femregion, Solution.t, not(round(t/Data.dt)==Data.AdaptivityStep));
           
        %% Project the solutions on the new femregion
        [Solution] = ProjectSolutionAssemblyFHN(Data, mesh, femregion, t, Solution);
    
        %% Reassemble the Matrices
        fprintf('\nMatrices computation ... \n');
        tic
        [Matrices] = IPDGMatrixAssemblyFHN(Data, mesh, femregion, Matrices);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')

        %% Right-hand side assembly
        fprintf('\nComputing forcing terms ... \n');
        [F_old] = ForcingTermAssemblyFHN(Data, mesh, femregion, t-Data.dt, Solution.t_old, F_old);
    
        [F_new] = ForcingTermAssemblyFHN(Data, mesh, femregion, t, Solution.t, F_new);
        
        %% Problem resolution with (possibly) nonlinear solver
        SolutionApp = PicardIterations(Data, mesh, femregion, Matrices, Solution.t_old, F_old, F_new, t);

        Solution.t.u_h = SolutionApp.u_h;
        Solution.t.w_h = SolutionApp.w_h;
        
        AdaptCount = AdaptCount + 1;
        
    end
end
