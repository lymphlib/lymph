%> @file  RunpConvergenceElastodynamics.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Convergence analysis for the Elastodynamics problem (polynomial refinements)
%>
%==========================================================================
%> @section classRunpConvergenceElastodynamics Class description
%==========================================================================
%> @brief          Sequence of run of MainEla.m
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

%% Simulation - Setup
run("../RunSetup.m")

%% Input Data - Boundary conditions - Forcing term
DataConvTestEla;

degree_vector = [1 2 3 4 5];

Errors.err_L2_d = [];
Errors.err_L2_v = [];
Errors.err_dG = [];
Errors.err_Energy = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq(1));
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfile,'P');
    end
    
    
    %% Main     
    [Error] = MainEla(Data,Setup);
    Errors.err_L2_v   = [Errors.err_L2_v,   Error.error_L2_v];
    Errors.err_L2_d   = [Errors.err_L2_d,   Error.error_L2_d];
    Errors.err_dG     = [Errors.err_dG,     Error.error_dG];
    Errors.err_Energy = [Errors.err_Energy, Error.error_Energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
semilogy(degree_vector,Errors.err_L2_d,'g')
hold on
semilogy(degree_vector,Errors.err_dG,'r')
semilogy(degree_vector,Errors.err_Energy,'b')
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
