%> @file  RunMainEla.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Run of MainEla for the solution of the elastodynamics problem
%>
%==========================================================================
%> @section classRunMainEla Class description
%==========================================================================
%> @brief            Run of MainEla
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
DataTestEla;
%DataTestPhysicsEla;

%% Mesh Generation
if Data.MeshFromFile
    % Load existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    % Create a new mesh
    Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
end

%% Main
[Error] = MainEla(Data,Setup);


