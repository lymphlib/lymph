%> @file  RunMainLaplacian.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Run of MainLaplacian for the solution of the Poisson problem
%>
%==========================================================================
%> @section classRunMainLaplacian Class description
%==========================================================================
%> @brief            Run of MainLaplacian
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================

%% Initial Simulation Setup
MyPhysicsPath = pwd;
run('../SimulationSetup.m');

%% Input Data - Boundary conditions - Forcing term
DataTestTanh;
%DataTestLap;

%% Mesh Generation

if Data.MeshFromFile
    % Load existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    % Create a new mesh
    Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
end

%% Main
[Error] = MainLaplacian(Data,Setup)
