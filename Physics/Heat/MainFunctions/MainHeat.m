%> @file  MainHeat.m
%> @author Mattia Corti
%> @date 7 July 2023
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
    fprintf('\n------------------------------------------------------------------\n')

    %% Load Region
    fprintf('\nLoading Region ... \n');

    mesh = load(Data.meshfile);

    % Compute mesh size
    Data.h = max(sqrt((mesh.region.BBox(:,1)-mesh.region.BBox(:,2)).^2  ...
        + (mesh.region.BBox(:,3)-mesh.region.BBox(:,4)).^2));

    fprintf(['Number of Polygonal Elements: ', num2str(mesh.region.ne)]);
    fprintf(['\nMesh size: ', num2str(Data.h)]);
    fprintf('\n\n------------------------------------------------------------------\n')


    % checking tags for elements
    for i = 1 : length(mesh.region.id)
        for j = 1 : length(Data.TagElLap)
            if mesh.region.id(i) == Data.TagElLap(j)
                mesh.region.tag(i,1) = 'L';
            end
        end
    end


    %% Create output directories
    if ( Setup.isSaveSolution || Setup.isSaveCSV ) && ~exist(Setup.OutFolder,'dir')
        mkdir(char(Setup.OutFolder));
    end
    if Setup.isSaveVTK && ~exist(Setup.OutFolderVTK,'dir')
        mkdir(char(Setup.OutFolderVTK));
    end


    %% Creation of the finite element space
    fprintf('\nMake femregion ... ')

    [femregion] = CreateDOF(Data, mesh.region);

    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')


    %% Plot polygonal mesh
    if Setup.isPlotMesh
        PlotPolymesh(mesh.neighbor, femregion);
    end

    %% Save VTK polygonal mesh
    if Setup.isSaveVTKMesh
        CreatePolygonalVTK(Data, Setup, mesh.region);
    end

    %% Matrix Assembly

    fprintf('\nMatrices computation ... \n');

    [Matrices] = MatrixHeat(Data, mesh.neighbor, femregion);

    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    %% Initial condition

    [u_old]  = EvaluateSolution(Data, femregion, Data.t0);

    % Projection for modal coordinates
    u_old  = Matrices.M_prj\u_old;


    %% Plot of initial condition

    [Gu] = GetSolutionQuadPoints(Data, femregion, mesh.neighbor, u_old, Data.t0);
    ScatterPlot(Gu, mesh.region, Data, 'u_0');

    %% Time-Loop Solver

    fprintf(['\nStarting Time: ', num2str(Data.t0)]);
    fprintf('\n------------------------------------------------------------------\n')

    [u_h, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, u_old);


    %% Computation of numerical errors

    if Setup.isError == 1

        fprintf('\nError computations...');
        fprintf('\n------------------------------------------------------------------\n')

        [Error] = ComputeErrors(Data, Matrices, femregion, u_h, Data.T);
        Error.err_Energy = sqrt(Error.err_u_L2.^2 + IntError);

    end

    Error.h = Data.h;


    %% Plot of final solution

    [Gu] = GetSolutionQuadPoints(Data, femregion, mesh.neighbor, u_h, Data.T);
    ScatterPlot(Gu, mesh.region, Data);


end
