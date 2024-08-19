%> @file  RunpConvergenceElastodynamics.m
%> @author The Lymph Team
%> @date 24 July 2024
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

Errors.err_u_L2   = [];
Errors.err_u_dG   = [];
Errors.err_v_L2   = [];
Errors.err_energy = [];
Errors.h          = [];

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
    Errors.err_u_L2   = [Errors.err_u_L2,   Error.err_u_L2];
    Errors.err_u_dG   = [Errors.err_u_dG,   Error.err_u_dG];
    Errors.err_v_L2   = [Errors.err_v_L2,   Error.err_v_L2];
    Errors.err_energy = [Errors.err_energy, Error.err_energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
semilogy(degree_vector,Errors.err_u_L2,'g')
hold on
semilogy(degree_vector,Errors.err_u_dG,'r')
semilogy(degree_vector,Errors.err_energy,'b')
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
