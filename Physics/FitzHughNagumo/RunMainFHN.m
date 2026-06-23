%> @file  RunMainFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 June 2026
%> @brief Run of MainFHN for the solution of the FHN problem
%>
%==========================================================================
%> @section classRunMainFHN Class description
%==========================================================================
%> @brief            Run of MainFHN
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
%DataTestCaseFHN
DataSpiral;

%% Mesh Generation
if Data.meshfromfile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.foldername,Data.meshfileseq);
else
    [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.foldername,Data.meshfileseq,'C');
    %[Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.foldername,Data.meshfileseq,'P');
end



%% Main 
[Error] = MainFHN(Data,Setup)
