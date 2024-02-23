%> @file  RunhConvergenceStokesPS.m
%> @author The Lymph Team
%> @date 16 Febraury 2024
%> @brief Convergence analysis for the unsteady Stokes problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergenceStokesPS Class description
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

Errors.err_dG = [];
Errors.err_Energy = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:4

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq(ii));
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq{ii},'P','laplacian');
    end
    
    
    %% Main     
    [Error] = MainStokesPS(Data,Setup);
    Errors.err_dG     = [Errors.err_dG,     Error.error_dG];
    Errors.err_Energy = [Errors.err_Energy, Error.error_Energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
Degree = Error.p;
figure(300);
loglog(Errors.h,Errors.err_Energy, Errors.h,Errors.h.^(Degree),'LineWidth',2);
legend('Error Energy norm', 'h^{p}')
