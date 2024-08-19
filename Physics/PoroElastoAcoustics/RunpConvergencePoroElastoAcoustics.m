%> @file  RunpConvergencePoroElastoAcoustics.m
%> @author The Lymph Team
%> @date 1 July 2024
%> @brief Convergence analysis for the coupled poro-elasto-acoustic problem (p refinements)
%>
%==========================================================================
%> @section classRunpConvergencePoroElastoAcoustics Class description
%==========================================================================
%> @brief          Sequence of run of RunpConvergencePoroElastoAcoustics.m
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
% DataConvTestAcu;
% DataConvTestEla;
% DataConvTestPoro;
DataConvTestPoroAcu;
% DataConvTestElaAcu
% DataConvTestPoroEla

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
        if (Data.NPhys == 1)
            Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','waves');
        elseif (Data.NPhys == 2)
            Data.meshfile = MakeMeshBidomainVert(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','waves');
        end
    end
    
    
    %% Main     
    [Error] = MainPoroAcuEla(Data,Setup);
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
