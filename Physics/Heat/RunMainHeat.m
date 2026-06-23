%> @file  RunMainHeat.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Run of MainHeat for the solution of the heat equation
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


%% Initial Simulation Setup
MyPhysicsPath = pwd;
run('../SimulationSetup.m');

%% Input Data - Boundary conditions - Forcing term
%DataTestCasepAdaptive;
DataTestCasePaper;

%% Mesh Generation

if Data.MeshFromFile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
    %Data.meshfile = MakeMeshDoubleCircle(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq);
end

%% Main 
[Error] = MainHeat(Data,Setup)
