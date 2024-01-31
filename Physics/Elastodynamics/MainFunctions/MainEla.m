%> @file  MainEla.m
%> @author Ilario Mazzieri
%> @date 16 September 2023
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

%% Check if output dir exists
if ~exist(Setup.OutFolder,'dir'); mkdir(Setup.OutFolder); end
if ~exist(Setup.OutFolderVTK,'dir'); mkdir(Setup.OutFolderVTK); end

fprintf('\nSolution of wave propagation problems \n');
fprintf('... here display some info about the simulation ... \n');
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
    
    %elastic element
    for j = 1 : length(Data.TagElEla)
        if mesh.region.id(i) == Data.TagElEla(j)
            mesh.region.tag(i,1) = 'E';
        end
    end
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


%%  Matrix Assembly
fprintf('\nMatrices computation ... \n');
tic
[Matrices] = MatWaves(Data, mesh.neighbor, femregion);
toc
fprintf('Done \n')
fprintf('\n------------------------------------------------------------------\n')


%% Right-hand side assembly

fprintf('\nComputing RHS ... \n');
tic
[F] = ForWaves(Data, mesh.neighbor, femregion);
toc
fprintf('Done\n')
fprintf('\n------------------------------------------------------------------\n')


%% Assembly of the ODE system A \ddot{x} + B\dot{x} + C x

fprintf('\nAssembly ODE system matrices ... \n');

[A, B, C] = AssembleODEMatrices(Matrices);

fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Evaluate initial conditions

fprintf('\nComputing initial conditions ... \n');

[Uold,Solutions] = GetInitalConditions(Data,femregion, Matrices);

%% Plot initial condition
if Data.PlotIniCond
    fprintf('\nPlot initial conditions ... \n');
    [Displacement, ~, Velocity] = GetSolutionQuadPoints(Data, femregion, Solutions, 0);
    ScatterPlot(Displacement.Ue, mesh.region, Data, Displacement.UeTag, 1); 
    ScatterPlot(Displacement.Ue, mesh.region, Data, Displacement.UeTag, 2); 
    ScatterPlot(Velocity.Ue, mesh.region, Data, Velocity.UeTag, 1); 
    ScatterPlot(Velocity.Ue, mesh.region, Data, Velocity.UeTag, 2); 
end


fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

if ~Setup.isError; clear Matrices; end

if strcmp(Data.timeint,'newmark')
    
    fprintf('\nNewmark time integration ... \n');
    tic    
    [Solutions] = NewmarkScheme(Setup, Data, femregion, mesh.region, A, B, C, F, Uold);
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    
else 
    
    fprintf('\nTime integration scheme undefined ... \n');  
   
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    stop;    
end

%% Plot solutions Matlab
if Setup.isPlotSolution
    fprintf('\nCompute Solution at quadrature points ... \n');
    
    [Displacement, ~, Velocity] = GetSolutionQuadPoints(Data, femregion, Solutions, Data.T);
    
    % Plot displacement, velocity and pressure
    fprintf('\nPlot Solution ... \n');
    ScatterPlot(Displacement.Ue, mesh.region, Data, Displacement.UeTag, 1); 
    ScatterPlot(Displacement.Ue, mesh.region, Data, Displacement.UeTag, 2); 
    ScatterPlot(Velocity.Ue, mesh.region, Data, Velocity.UeTag, 1); 
    ScatterPlot(Velocity.Ue, mesh.region, Data, Velocity.UeTag, 2); 
    
    %PlotPressure(Pressure,20)
    
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
end



%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');
    
    [Error] = ComputeErrors(Data, femregion, Matrices, Solutions);
    
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

fprintf('\nBye \n');

end

