%> @file  RunMainFKPP.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Run of MainFKPP for the solution of the Fisher-KPP problem
%>
%==========================================================================
%> @section classRunMainFKPP Class description
%==========================================================================
%> @brief            Run of MainFKPP
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
DataWavesFKPP;

%% Mesh Generation
if Data.meshfromfile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.foldername,Data.meshfileseq);
else
    [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.foldername,Data.meshfileseq,'P');
end

%% Main 
[Error] = MainFKPP(Data,Setup)
