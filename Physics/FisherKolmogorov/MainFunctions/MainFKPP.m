%> @file  MainFKPP.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Solution of the FisherKPP problem with PolydG
%>
%==========================================================================
%> @section classMainFKPP Class description
%==========================================================================
%> @brief            Solution of the Fisher KPP problem with PolydG
%>
%> @param Data       Struct with problem's data
%> @param Setup      Simulation's setup
%>
%> @retval Error     Struct for convergence test
%>
%==========================================================================

function [Error] = MainFKPP(Data, Setup)

    fprintf('\nSolution of the Fisher-KPP problem \n');
    fprintf('... here display some info about the simulation ... \n');
    fprintf('\n------------------------------------------------------------------\n')

    [mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data);

    %% Matrix Assembly

    fprintf('\nMatrices computation ... \n');
    
    [Matrices] = IPDGMatrixAssemblyFKPP(Data, mesh, femregion, []);
    
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    
    %% Initial conditions evaluation
    Solution.c_h = ComputeModalSolutionFKPP(Data, mesh, femregion, Data.t0);
    Solution.c_h = Matrices.M_prj\Solution.c_h;

    Solution.c_old = ComputeModalSolutionFKPP(Data, mesh, femregion, Data.t0-Data.dt);
    Solution.c_old = Matrices.M_prj\Solution.c_old;

    %% Compute the indicator for p-adaptivity
    if Data.Adaptivity
        fprintf('\nComputing adaptive indicator ... ');
        tic
        [Solution.Indicator] = ComputeAdaptiveIndicatorFKPP(Data, mesh, femregion, Solution, Data.t0);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end

    %% Postprocessing of the initial condition
    if Data.PlotIniCond
        fprintf('\nPlot initial conditions ... \n');
        PostProcessSolution(Setup, Data, mesh, femregion, 0, Solution, Data.t0);
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end

    %% Time-Loop Solver
    fprintf(['\nStarting Time: ', num2str(Data.t0)]);
    fprintf('\n------------------------------------------------------------------\n')
    
    [Solution, femregion, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, Solution.c_h);
    
    %% Computation of numerical errors

    if Setup.isError == 1
    
        fprintf('\nError computations...');
        fprintf('\n------------------------------------------------------------------\n')
        [Error] = ComputeErrorsFKPP(Data, mesh, femregion, Solution.c_h, Data.T);
        Error.err_Energy = sqrt(Error.L2.^2 + IntError);

        Error.nel = femregion.nel;
        Error.h   = Data.h;
        Error.p   = Data.degree;
    end

    Error.h = Data.h;

end
