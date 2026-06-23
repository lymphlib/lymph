%> @file  MainPoroAcuEla.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 5 June 2026
%> @brief Solution of the coupled poro-elasto-acoustic problem with PolydG
%>
%==========================================================================
%> @section classMainPoroAcuEla Class description
%==========================================================================
%> @brief            Solution of the poro-elasto-acoustic problem with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation setup
%
%> @retval Error     Struct with computed errors
%>
%==========================================================================

function [Error] = MainPoroAcuEla(Data,Setup)

fprintf('\nSolution of poro-elasto-acostic propagation problems \n');
fprintf('... here display some info about the simulation ... \n');
fprintf('\n------------------------------------------------------------------\n')

%% Load Region
[mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data);

%%  Matrix Assembly
fprintf('\nMatrices computation ... \n');
tic
[Matrices] = MatrixAssemblyPoroElastoAcoustics(Data, mesh, femregion);
toc
fprintf('Done \n')
fprintf('\n------------------------------------------------------------------\n')


%% Right-hand side assembly

fprintf('\nComputing RHS ... \n');
tic
[F] = ForcingTermAssemblyPoroElastoAcoustics(Data, mesh, femregion, Data.t0);
toc
fprintf('Done\n')
fprintf('\n------------------------------------------------------------------\n')

%% Assembly of the ODE system A \ddot{x} + B\dot{x} + C x

fprintf('\nAssembly ODE system matrices ... \n');
[A, B, C] = AssembleODEMatricesPoroAcuEla(Data,Matrices);
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Evaluate initial conditions

fprintf('\nComputing initial conditions ... \n');

[Uold, Solutions] = GetInitalConditionsPoroAcuEla(Data, mesh, femregion, Matrices);

%% Plot initial condition
if Data.PlotIniCond
    fprintf('\nPlot initial conditions ... \n');
    PostProcessSolution(Setup, Data, mesh, femregion, 0, Solutions, Data.t0);
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
end
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

if ~Setup.isError; clear Matrices; end

if strcmp(Data.timeint,'newmark')

    fprintf('\nNewmark time integration ... \n');
    tic
    [Solutions] = NewmarkSchemePoroAcuEla(Setup, Data, femregion, mesh, A, B, C, F, Uold);
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

else

    fprintf('\nTime integration scheme undefined ... \n');

    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    stop;
end

%% Post-processing final time

if Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution
    PostProcessSolution(Setup, Data, mesh, femregion, floor(Data.T/Data.dt)+1, Solutions, Data.T);
end

%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');

    [Error] = ComputeErrorsPoroAcuEla(Data, mesh, femregion, Solutions, Data.T);
    Error.h = Data.h;
    
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

fprintf('\nBye \n');

