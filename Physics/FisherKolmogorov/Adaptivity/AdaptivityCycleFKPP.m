%> @file  AdaptivityCycleFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Implementation of adaptivity cycle for the FKPP equation.
%>
%==========================================================================
%> @section classAdaptivityCycleFKPP Class description
%==========================================================================
%> @brief Implementation of adaptivity cycle for the FKPP equation.
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

function [Data, Solution, femregion, Matrices, F_new] = AdaptivityCycleFKPP(Data, mesh, femregion, Solution, Matrices, F_old, F_new, t)

    AdaptCount = 0;
    
    while AdaptCount < Data.AdaptIts
    
        fprintf('\nComputing adaptive indicator ... ');
        tic
        %% Indicator computation
        [Solution.Indicator] = ComputeAdaptiveIndicatorFKPP(Data, mesh, femregion, Solution, t);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    
        %% Update the polynomial degrees
        [Data, femregion] = ComputeUpdatedDegrees(Data, mesh, femregion, Solution, not(round(t/Data.dt)==Data.AdaptivityStep));
           
        %% Project the solutions on the new femregion
        [Solution] = ProjectSolutionAssemblyFKPP(Data, mesh, femregion, t, Solution);
    
        %% Reassemble the Matrices
        fprintf('\nMatrices computation ... \n');
        tic
        [Matrices] = IPDGMatrixAssemblyFKPP(Data, mesh, femregion, Matrices);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')

        %% Right-hand side assembly
        fprintf('\nComputing forcing terms ... \n');
        [F_old] = ForcingTermAssemblyFKPP(Data, mesh, femregion, t-Data.dt, F_old);
    
        [F_new] = ForcingTermAssemblyFKPP(Data, mesh, femregion, t, F_new);
        
        %% Problem resolution with (possibly) nonlinear solver
        Solution.c_h = PicardIterations(Data, femregion, Matrices, Solution, F_old, F_new);

        AdaptCount = AdaptCount + 1;
    end
end
