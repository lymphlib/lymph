%> @file  MainStokesPS.m
%> @author Ilario Mazzieri
%> @date 13 February 2024
%> @brief Solution of the unseteady Stokes equation with PolydG
%>
%==========================================================================
%> @section classMainStokesPS Class description
%==========================================================================
%> @brief            Solution of the Stokes eq. with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation setup
%
%> @retval Error     Struct with computed errors
%>
%==========================================================================
function [Error] = MainStokesPS(Data,Setup)

%% Check if output dir exists
if ~exist(Setup.OutFolder,'dir'); mkdir(Setup.OutFolder); end
if ~exist(Setup.OutFolderVTK,'dir'); mkdir(Setup.OutFolderVTK); end

fprintf('\nSolution of Stokes problem in pseudo-stress formulation \n');
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

    % fluid element
    for j = 1 : length(Data.TagElFluid)
        if mesh.region.id(i) == Data.TagElFluid(j)
            mesh.region.tag(i,1) = 'F';
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
    [~] = PlotPolymesh(mesh.neighbor, femregion);
end

%% Save VTK polygonal mesh
if Setup.isSaveVTKMesh
    CreatePolygonalVTK(Data, Setup, mesh.region);
end

%%  Matrix Assembly
fprintf('\nMatrices computation ... \n');

[Matrices] = MatStokesPS(Data, mesh.neighbor, femregion);

fprintf('Done \n')
fprintf('\n------------------------------------------------------------------\n')


%% Right-hand side assembly

fprintf('\nComputing RHS ... \n');

[F] = ForStokesPS(Data, mesh.neighbor, femregion);

fprintf('Done\n')
fprintf('\n------------------------------------------------------------------\n')


%% Assembly of the ODE system M \dot{x} + K x = F

fprintf('\nAssembly ODE system matrices ... \n');

M = Matrices.Fluid.M_F;
K = Matrices.Fluid.A_F;

fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Evaluate initial conditions

fprintf('\nComputing initial conditions ... \n');

[Uold,Solutions] = GetInitalConditions(Data,femregion, Matrices);

% [Displacement, Pressure, Velocity] = GetSolutionQuadPoints(Data, femregion, Solutions, 0);
% PlotSolution(Data,Displacement,10);
% PlotSolution(Data,Velocity,30)

%% Plot initial condition
if Data.PlotIniCond
    fprintf('\nPlot initial conditions ... \n');
    [Disp, ~] = GetSolutionQuadPointsPS(Data, femregion, Solutions, 0);
    ScatterPlot(Disp.Sigma(:,[1 2 3 7]), mesh.region, Data, [Disp.SigmaTag,'11'], 0);
    ScatterPlot(Disp.Sigma(:,[1 2 4 8]), mesh.region, Data, [Disp.SigmaTag,'12'], 0);
    ScatterPlot(Disp.Sigma(:,[1 2 5 9]), mesh.region, Data, [Disp.SigmaTag,'21'], 0);
    ScatterPlot(Disp.Sigma(:,[1 2 6 10]), mesh.region, Data, [Disp.SigmaTag,'22'], 0);
end


fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

if ~Setup.isError; clear Matrices; end

if strcmp(Data.timeint,'CN')

    fprintf('\nCrank-Nicolson time integration ... \n');

    [Solutions] = CNScheme(Setup, Data, mesh, femregion, M, K, F, Uold);

    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

else

    fprintf('\nTime integration scheme undefined ... \n');

    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    stop;
end


%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');

    [Error] = ComputeErrorsStokesPS(Data, femregion, Matrices, Solutions);

    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

%% Compute Pressure and Velocity fields
if Data.ComputeVelAndPres
    fprintf('\nPost Processing Solution to compute Velocity and Pressure... \n');

    [Vel,Pres] = ComputeVelocityAndPressure(Data, femregion, Solutions);
    csvwrite(fullfile(Setup.OutFolder,'Vel.csv'),Vel);
    csvwrite(fullfile(Setup.OutFolder,'Pres.csv'),Pres);


    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
end


fprintf('\nBye \n');

end

