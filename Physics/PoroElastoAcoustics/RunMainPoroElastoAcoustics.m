%> @file  RunMainPoroElastoAcoustics.m
%> @author The Lymph Team
%> @date 15 June 2024
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
% DataTestAcu;             %Convergence test acoustic - Dirichlet
% DataTestEla;             %Convergence test elastic - Dirichlet
% DataTestPoro;            %Convergence test poroelastic - Dirichlet
% DataTestPoroAcu;         % Convergence test poroacoustic - paper SISC
% DataTestPoroAcuHor;         % Convergence test poroacoustic - paper SISC
% DataTestElaAcu;          % Convergence test elastoacoustic - paper CMAME
% DataTestPoroEla;         % Convergence test poro-elastic - paper IMAJNA (submitted)
DataTestPoroAcuEla;        % Convergence test poro-elastic-acoustic

% DataTestPhysicsPoroAcuEla;


%% Mesh Generation
if Data.MeshFromFile
    % Load existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    % Create a new mesh
    if (Data.NPhys == 1)
        Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','waves');
    elseif (Data.NPhys == 2)
        % Data.meshfile = MakeMeshBidomainVert(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','waves');
        Data.meshfile = MakeMeshBidomainHor(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','waves');
    elseif (Data.NPhys == 3)
        Data.meshfile = MakeMeshQuadridomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','waves');
    end
end

%% Main
[Error] = MainPoroAcuEla(Data,Setup);


