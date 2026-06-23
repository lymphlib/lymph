%> @file  MainFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 June 2026
%> @brief Solution of the FHN problem with PolydG
%>
%==========================================================================
%> @section classMainFHN Class description
%==========================================================================
%> @brief            Solution of the FHN problem with PolydG
%>
%> @param Data       Struct with problem's data
%> @param Setup      Simulation's setup
%>
%> @retval Error     Struct for convergence test
%>
%==========================================================================

function [Error] = MainFHN(Data, Setup)

    fprintf('\nSolution of the FHN problem \n');
    fprintf('... here display some info about the simulation ... \n');
    fprintf('\n------------------------------------------------------------------\n')

    [mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data);

    %% Matrix Assembly

    fprintf('\nMatrices computation ... \n');
    
    [Matrices] = IPDGMatrixAssemblyFHN(Data, mesh, femregion, []);
    
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    
    %% Initial conditions evaluation
    [Solution.t] = ComputeModalSolutionFHN(Data, mesh, femregion, Data.t0);
    Solution.t.u_h = Matrices.M_prj\Solution.t.u_h;
    Solution.t.w_h = Matrices.M_prj\Solution.t.w_h;
    
    [Solution.t_old] = ComputeModalSolutionFHN(Data, mesh, femregion, Data.t0-Data.dt);
    Solution.t_old.u_h = Matrices.M_prj\Solution.t_old.u_h;
    Solution.t_old.w_h = Matrices.M_prj\Solution.t_old.w_h;

    %% Compute the indicator for p-adaptivity
    if Data.Adaptivity
        fprintf('\nComputing adaptive indicator ... ');
        tic
        [Solution.t.Indicator] = ComputeAdaptiveIndicatorFHN(Data, mesh, femregion, Solution, Data.t0);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end

    %% Postprocessing of the initial condition
    if Data.PlotIniCond
        fprintf('\nPlot initial conditions ... \n');
        PostProcessSolution(Setup, Data, mesh, femregion, 0, Solution.t, Data.t0);
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end

    %% Time-Loop Solver
    fprintf(['\nStarting Time: ', num2str(Data.t0)]);
    fprintf('\n------------------------------------------------------------------\n')
    
    [Solution, femregion, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, Solution);
    
    %% Computation of numerical errors

    if Setup.isError
        fprintf('\nComputing Errors ... \n');
        [Error] = ComputeErrorsFHN(Data, mesh, femregion, Solution.t.u_h, Data.T);

        Error.nel = femregion.nel;
        Error.h   = Data.h;
        Error.p   = Data.degree;
        Error.NDoF= femregion.ndof;

        Solution.Cells_L2 = Error.Cells_L2;
        Solution.Cells_dG = Error.Cells_dG;

        fprintf('Done\n')
        fprintf('\n------------------------------------------------------------------\n')
    else
        Error = [];
    end

end
