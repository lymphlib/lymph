%> @file  RunMainHeat.m
%> @author The Lymph Team
%> @date 9 October 2023
%> @brief Run of MainLaplacian for the solution of the heat equation
%>
%==========================================================================
%> @section classRunMainHeat Class description
%==========================================================================
%> @brief            Run of MainHeat
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
DataTestCasePaper;

%% Mesh Generation

if Data.MeshFromFile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    Data.meshfile = MakeMeshDoubleCircle(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'laplacian');
end



%% Main 

[Error] = MainHeat(Data,Setup)
