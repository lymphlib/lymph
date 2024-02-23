%> @file  RunhConvergenceStokesPS.m
%> @author The Lymph Team
%> @date 16 Febraury 2024
%> @brief Convergence analysis for the unsteady Stokes problem (mesh refinements)
%>
%==========================================================================
%> @section classRunpConvergenceStokesPS Class description
%==========================================================================
%> @brief          Sequence of run of MainStokesPS.m
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================

%% Import lymph and paths of folders related to this problem
run("../ImportLymphPaths.m")
MyPhysicsPath = pwd;
addpath(genpath(fullfile(MyPhysicsPath,'Assembly')));
addpath(genpath(fullfile(MyPhysicsPath,'Error')));
addpath(genpath(fullfile(MyPhysicsPath,'InputData')));
addpath(genpath(fullfile(MyPhysicsPath,'InputMesh')));
addpath(genpath(fullfile(MyPhysicsPath,'MainFunctions')));
addpath(genpath(fullfile(MyPhysicsPath,'PostProcessing')));
addpath(genpath(fullfile(MyPhysicsPath,'TimeIntegration')));
addpath(genpath(fullfile(MyPhysicsPath,'Utilities')));

%% Simulation - Setup
run("../RunSetup.m")

%% Input Data - Boundary conditions - Forcing term
DataConvTestFluid;

degree_vector = [1 2 3 4 5];
Errors.err_dG = [];
Errors.err_Energy = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,'UnitSquareMixed_50_el.mat');
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfile,'P');
    end
    
    
    %% Main     
    [Error] = MainStokesPS(Data,Setup);
    Errors.err_dG     = [Errors.err_dG,     Error.error_dG];
    Errors.err_Energy = [Errors.err_Energy, Error.error_Energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
semilogy(degree_vector,Errors.err_Energy,'b', 'LineWidth',2)
legend("Error energy-norm","Interpreter","latex")
grid on
