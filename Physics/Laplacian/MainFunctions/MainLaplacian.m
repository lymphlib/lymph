%> @file  MainLaplacian.m
%> @author The Lymph Team
%> @date 26 July 2024
%> @brief Solution of the Poisson problem with PolydG
%>
%==========================================================================
%> @section classMainLaplacian Class description
%==========================================================================
%> @brief            Solution of the Poisson problem with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation's setup
%>
%> @retval Error     Struct for convergence test
%>
%==========================================================================

function [Error] = MainLaplacian(Data,Setup)

fprintf('\nSolution of the Poisson problem \n');
fprintf('... here display some info about the simulation ... \n');
fprintf('\n------------------------------------------------------------------\n')

%% Mesh setup
[mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data, {Data.TagElLap}, {'L'});

%% Matrix Assembly
fprintf('\nMatrices computation ... \n');
tic
switch Data.quadrature
    case "QF"
        [Matrices] = MatrixLaplacianQF(Data, mesh.neighbor, femregion);
    case "ST"
        [Matrices] = MatrixLaplacianST(Data, mesh.neighbor, femregion);
end
toc
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Right-hand side assembly

fprintf('\nComputing RHS ... \n');
tic
[F] = ForcingLaplacian(Data, mesh.neighbor, femregion);
toc
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Solving the linear system
fprintf('\nSolving linear system ... ');
tic
U = Matrices.A \ F;
toc
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Postprocess solution
if Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution
    PostProcessSolution(Setup, Data, mesh, femregion, 0, U);
end

%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');
    [Error] = ComputeErrors(Data, femregion, Matrices, U);
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

fprintf('\nBye \n');


