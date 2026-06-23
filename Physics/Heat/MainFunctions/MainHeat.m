%> @file  MainHeat.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Solution of the heat equation with PolydG
%>
%==========================================================================
%> @section classHeatMainHeat Class description
%==========================================================================
%> @brief            Solution of the heat equation with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation's setup
%>
%> @retval Error     Struct for convergence test
%>
%==========================================================================

function [Error] = MainHeat(Data, Setup)

    fprintf('\nSolution of the heat equation \n');
    fprintf('... here display some info about the simulation ... \n');
    fprintf('\n------------------------------------------------------------------\n')

    [mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data);

    %% Matrix Assembly

    fprintf('\nMatrices computation ... \n');
    tic
    [Matrices] = MatrixAssemblyHeat(Data, mesh, femregion, []);
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    %% Initial condition
    fprintf('\nComputing initial conditions ... \n');
    [u_h] = ComputeModalSolutionHeat(Data, mesh, femregion, Data.t0);
    Solution.u_h = Matrices.Mprj\u_h;
    [u_old] = ComputeModalSolutionHeat(Data, mesh, femregion, Data.t0-Data.dt);
    Solution.u_old = Matrices.Mprj\u_old;
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    %% Compute the indicator for p-adaptivity
    if Data.Adaptivity
        fprintf('\nComputing adaptive indicator ... ');
        tic
        [Solution.Indicator] = ComputeAdaptiveIndicatorHeat(Data, mesh, femregion, Solution, Data.t0);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end

    %% Plot of initial condition
    if Data.PlotIniCond
        fprintf('\nPlot initial conditions ... \n');
        PostProcessSolution(Setup, Data, mesh, femregion, 0, Solution, Data.t0);
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    %% Time-Loop Solver

    fprintf(['\nStarting Time: ', num2str(Data.t0)]);
    fprintf('\n------------------------------------------------------------------\n')

    [Solution, femregion, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, Solution.u_h);

    %% Computation of numerical errors
    if Setup.isError == 1

        fprintf('\nError computations...');
        fprintf('\n------------------------------------------------------------------\n')
        [Error] = ComputeErrorsHeat(Data, mesh, femregion, Solution.u_h, Data.T);
        Error.err_Energy = sqrt(Error.L2.^2 + IntError);

        Error.nel = femregion.nel;
        Error.h   = Data.h;
        Error.p   = Data.degree;
    else
        Error = [];
    end

    fprintf('\nBye \n');

end
