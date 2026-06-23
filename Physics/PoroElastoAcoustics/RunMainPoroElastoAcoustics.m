%> @file  RunMainPoroElastoAcoustics.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Run of MainPoroElastoAcoustics for the solution of the coupled
%>   poroelastic-elstic-acoustic problem
%>
%==========================================================================
%> @section classRunMainPoroElaAcu Class description
%==========================================================================
%> @brief            Run of MainPoroElaAcu
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================

%% Initial Simulation Setup
MyPhysicsPath = pwd;
run('../SimulationSetup.m');

% Reuse of the assembly routines of elastodynamics (excluded faces assembly) to avoid code duplication
addpath(genpath('..\Elastodynamics\Assembly\'));

%% Input Data - Boundary conditions - Forcing term
% DataTestAcu;             %Convergence test acoustic - Dirichlet
% DataTestEla;             %Convergence test elastic - Dirichlet
% DataTestPoro;            %Convergence test poroelastic - Dirichlet
% DataTestPoroAcu;         % Convergence test poroacoustic - paper SISC
% DataTestPoroAcuHor;      % Convergence test poroacoustic - paper SISC
% DataTestElaAcu;          % Convergence test elastoacoustic - paper CMAME
% DataTestPoroEla;         % Convergence test poro-elastic - paper IMAJNA (submitted)
% DataTestPoroAcuEla;      % Convergence test poro-elastic-acoustic

%% Mesh Generation
if Data.MeshFromFile
    % Load existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    % Create a new mesh
    if (Data.NPhys == 1)
        Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
    elseif (Data.NPhys == 2)
        % Data.meshfile = MakeMeshBidomainVert(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
        Data.meshfile = MakeMeshBidomainHor(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
    elseif (Data.NPhys == 3)
        Data.meshfile = MakeMeshQuadridomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
    end
end

%% Main
[Error] = MainPoroAcuEla(Data,Setup);


