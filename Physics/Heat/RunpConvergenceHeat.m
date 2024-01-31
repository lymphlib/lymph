%> @file  RunpConvergenceHeat.m
%> @author The Lymph Team
%> @date 9 October 2023
%> @brief Convergence analysis for the heat equation (polynomial refinements)
%>
%==========================================================================
%> @section classRunpConvergenceHeat Class description
%==========================================================================
%> @brief          Sequence of run of MainHeat.m
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================


%% Import lymph and add path related to this physics.
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
DatapConvergenceTest;

degree_vector = Data.degree;

Errors.err_u_L2 = [];
Errors.err_u_dG = [];
Errors.err_Energy = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq);
    else
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfile,'P');
    end
    
    
    
    %% Main 
    
    [Error] = MainHeat(Data,Setup);
    Errors.err_u_L2 = [Errors.err_u_L2, Error.err_u_L2];
    Errors.err_u_dG = [Errors.err_u_dG, Error.err_u_dG];
    Errors.err_Energy = [Errors.err_Energy, Error.err_Energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
semilogy(degree_vector,Errors.err_u_L2,'g')
hold on
semilogy(degree_vector,Errors.err_u_dG,'r')
semilogy(degree_vector,Errors.err_Energy,'b')
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
