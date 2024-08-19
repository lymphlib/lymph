%> @file  MainHeat.m
%> @author Mattia Corti
%> @date 26 July 2024
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

    [mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data, {Data.TagElLap}, {'L'});

    %% Matrix Assembly

    fprintf('\nMatrices computation ... \n');
    tic
    switch Data.quadrature
        case "QF"
            [Matrices] = MatrixHeatQF(Data, mesh.neighbor, femregion);
        case "ST"
            [Matrices] = MatrixHeatST(Data, mesh.neighbor, femregion);
    end
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    %% Initial condition
    fprintf('\nComputing initial conditions ... \n');
    [u_old] = GetInitalConditions(Data, femregion, Matrices);
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')


    %% Plot of initial condition

    if Data.PlotIniCond
        fprintf('\nPlot initial conditions ... \n');
        PostProcessSolution(Setup, Data, mesh, femregion, 0, u_old, Data.t0);
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')
    end
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    %% Time-Loop Solver

    fprintf(['\nStarting Time: ', num2str(Data.t0)]);
    fprintf('\n------------------------------------------------------------------\n')

    [u_h, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, u_old);


    %% Post-processing final time

    if Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution
        PostProcessSolution(Setup, Data, mesh, femregion, floor(Data.T/Data.dt)+1, u_h, Data.T);
    end

    %% Computation of numerical errors

    if Setup.isError == 1

        fprintf('\nError computations...');
        fprintf('\n------------------------------------------------------------------\n')
        [Error] = ComputeErrors(Data, Matrices, femregion, u_h, Data.T);
        Error.err_Energy = sqrt(Error.err_u_L2.^2 + IntError);
    else
        Error = [];
    end

    fprintf('\nBye \n');

end
