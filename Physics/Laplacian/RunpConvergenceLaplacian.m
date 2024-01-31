%> @file  RunpConvergenceLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Convergence analysis for the Poisson problem (polynomial refinement)
%>
%==========================================================================
%> @section classRunpConvergenceLaplacian Class description
%==========================================================================
%> @brief          Sequence of run of MainLaplacian.m
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
addpath(genpath(fullfile(MyPhysicsPath,'InputData')));
addpath(genpath(fullfile(MyPhysicsPath,'MainFunctions')));
addpath(genpath(fullfile(MyPhysicsPath,'Error')));
addpath(genpath(fullfile(MyPhysicsPath,'PostProcessing')));

%% Simulation - Setup
run("../RunSetup.m")

%% Input Data - Boundary conditions - Forcing term
DataConvTestLap;

degree_vector = [1 2 3 4 5 6];

Errors.err_L2 = [];
Errors.err_dG = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq{1});
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N{1},Data.domain,Data.FolderName,Data.meshfileseq{1},'P','laplacian');
    end
    
    
    %% Main     
    [Error] = MainLaplacian(Data,Setup);
    Errors.err_L2   = [Errors.err_L2, Error.L2];
    Errors.err_dG   = [Errors.err_dG, Error.dG];
    Errors.h        = [Errors.h, Error.h];

end

%% Plot of the errors
figure
semilogy(degree_vector,Errors.err_L2,'g','Linewidth',2)
hold on
semilogy(degree_vector,Errors.err_dG,'r','Linewidth',2)
legend("Error $L^2$-norm", "Error DG-norm", "Interpreter","latex")
xlabel('$\ell$',"Interpreter","latex")
grid on
