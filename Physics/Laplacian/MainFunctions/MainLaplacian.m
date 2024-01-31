%> @file  MainLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
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

%% Check if output dir exists
if ~exist(Setup.OutFolder,'dir'); mkdir(Setup.OutFolder); end
if ~exist(Setup.OutFolderVTK,'dir'); mkdir(Setup.OutFolderVTK); end


fprintf('\nSolution of the Poisson problem \n');
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
    for j = 1 : length(Data.TagElLap)
        if mesh.region.id(i) == Data.TagElLap(j)
            mesh.region.tag(i,1) = 'L';
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


%% Matrix Assembly
fprintf('\nMatrices computation ... \n');
tic
[Matrices] = MatrixLaplacian(Data, mesh.neighbor, femregion);
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
Solutions.U = Matrices.A \ F;
toc
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')


%% Computing solution at quadrature points
fprintf('\nComputing solution at quadrature nodes ... \n');
tic
[UPlot] = GetSolutionQuadPoints(Data, femregion, mesh.neighbor, Solutions);
toc
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')


%% Save .mat file
if Setup.isSaveSolution
    FileName = fullfile(Setup.OutFolder,[Data.name,'.mat']);
    save(FileName,'Data','femregion','Solutions');
end

%% Save solutions paraview
if Setup.isSaveCSV
    CurrentFolder = pwd;
    filenameout = fullfile(Setup.OutFolder,['Solution', '.csv']);

    if Data.PlotExact
        title = {'x','y','u_h','u_ex'};
    else
        title = {'x','y','u_h'};
    end
    UPlot = array2table(UPlot);
    UPlot.Properties.VariableNames(1:length(title)) = title;
    writetable(UPlot,filenameout);
end
if Setup.isSaveVTK
    PlotSolutionParaview(Setup, Data, UPlot, 1,'uh');
end

%% Plot solution matlab
if Setup.isPlotSolution
    ScatterPlot(UPlot, mesh.region, Data, 'u');
end


%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');
    tic 
    [Error] = ComputeErrors(Data, femregion, Matrices, Solutions);
    toc
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

fprintf('\nBye \n');


