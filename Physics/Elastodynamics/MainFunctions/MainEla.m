%> @file  MainEla.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti
%> @date 23 July 2024
%> @brief Solution of the elastodynamics equation with PolydG
%>
%==========================================================================
%> @section classMainEla Class description
%==========================================================================
%> @brief            Solution of the elastodynamics eq. with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation setup
%
%> @retval Error     Struct with computed errors 
%>
%==========================================================================

function [Error] = MainEla(Data,Setup)

fprintf('\nSolution of wave propagation problems \n');
fprintf('... here display some info about the simulation ... \n');
fprintf('\n------------------------------------------------------------------\n')

[mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data, {Data.TagElEla}, {'E'});

%%  Matrix Assembly

fprintf('\nMatrices computation ... \n');
tic
switch Data.quadrature
    case "QF"
        [Matrices] = MatElaQF(Data, mesh.neighbor, femregion);
    case "ST"
        [Matrices] = MatElaST(Data, mesh.neighbor, femregion);
end
toc
fprintf('Done \n')
fprintf('\n------------------------------------------------------------------\n')

%% Assembly of the ODE system A \ddot{x} + B\dot{x} + C x

fprintf('\nAssembly ODE system matrices ... \n');
[A, B, C] = AssembleODEMatrices(Matrices);
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Evaluate initial conditions

fprintf('\nComputing initial conditions ... \n');
U_old = GetInitalConditions(Data,femregion, Matrices);
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Plot initial condition

if Data.PlotIniCond
    fprintf('\nPlot initial conditions ... \n');
    PostProcessSolution(Setup, Data, mesh, femregion, 0, U_old, Data.t0);
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
end
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

if ~Setup.isError; clear Matrices; end

%% Time-integration loop

if strcmp(Data.timeint,'newmark')
    fprintf('\nNewmark time integration ... \n');
    tic    
    U_h = NewmarkScheme(Setup, Data, femregion, mesh, A, B, C, U_old);
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
    PostProcessSolution(Setup, Data, mesh, femregion, floor(Data.T/Data.dt)+1, U_h, Data.T);
end

%% Compute Errors

if Setup.isError    
    fprintf('\nComputing Errors ... \n');
    [Error] = ComputeErrors(Data, femregion, Matrices, U_h);
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else   
    Error = [];
end

fprintf('\nBye \n');

end