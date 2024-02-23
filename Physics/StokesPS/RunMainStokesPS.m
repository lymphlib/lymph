%> @file  RunMainStokesPS.m
%> @author The Lymph Team
%> @date 13 Febraury 2024
%> @brief Run of MainStokesPS for the solution of the unsteady Stokes
%problem in pseudo-stress formulation
%>
%==========================================================================
%> @section classRunMainStokesPS Class description
%==========================================================================
%> @brief            Run of RunMainStokesPS
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
DataTestFluid;
% DataTestVelRecFluid;
% DataTestFACFluid;

%% Mesh Generation
if Data.MeshFromFile
    % Load existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    % Create a new mesh
    Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','stokes');
end

[Error] = MainStokesPS(Data,Setup);


% ------------------------- Mesh Generation -------------------------------
% ------------------------ Define a new mesh ------------------------------
% Rectangular domain
% N = 200;
% DomainLimits = [0 1 0 1];
% FolderName = 'InputMesh';
% FileName = 'UnitSquareVelRec';
% [Data.meshfile] = MakeMeshMonodomain(Data,N,DomainLimits,FolderName,FileName,'P','stokes');

% %Rectangular domain with circular inclusion
% N = 2000;
% DomainLimits = [-1 4 -1 1];
% Circle.Radius = 0.2;
% Circle.Center = [0 0];
% FolderName = 'InputMesh/';
% FileName = 'FlowAroundCilinder';
% [Data.meshfile] = MakeMeshFlowAroundCilinder(Data,N,DomainLimits,Circle,FolderName,FileName);

% ------------------- Load an existing mesh -------------------------------
% Data.meshfile = fullfile('InputMesh','UnitSquareVelRec_200_el.mat'); % fluid
% Data.meshfile = fullfile('InputMesh','FlowAroundCilinder_2000_el.mat'); %fluid
% Data.meshfile = fullfile('InputMesh','UnitSquareMixed_100_el.mat'); % fluid




