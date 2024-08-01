%> @file  MainExample.m
%> @author Ilario Mazzieri
%> @date 8 March 2023
%> @brief Simple example of mesh generation

clear; clc; close all;
addpath(genpath('./Polymesher'))


%% 
N = 100;
DomainLimits = [-1 4 -1 1];
Circle.Radius = 0.2;
Circle.Center = [0 0];
FolderName = 'MeshTest/';
FileName = 'FlowAroundCilinder';
Data = struct;
Data.TagBcFluid = [2 3 4 5];
Data.LabBcFluid = 'DDDD';
[FileNameOutCircle] = MakeMeshFlowAroundCilinder(Data,N,DomainLimits,Circle,FolderName,FileName);



%%
Data = struct;

N = 100;
DomainLimits = [0 1 0 1];
FolderName = 'MeshTest/';
FileName = 'MonoDomainTest';
[FileNameOutMono] = MakeMeshMonodomain(Data,N,DomainLimits,FolderName,FileName,'P');

N = [10 10];
DomainLimits = [0 1 0 1];
FolderName = 'MeshTest/';
FileName = 'Quad01_quad_ref0';
[FileNameOutMonoC] = MakeMeshMonodomain(Data,N,DomainLimits,FolderName,FileName,'C');

%%
Data = struct;

N = 400;
DomainLimits = [-2400 2400 0 4800];
FolderName = 'MeshTest/';
FileName = 'BiDomainTestVer';
[FileNameOutBiVer] = MakeMeshBidomainVert(Data,N,DomainLimits,FolderName,FileName,'P');

N = [20 20];
DomainLimits = [-2400 2400 0 4800];
FolderName = 'MeshTest/';
FileName = 'BiDomainTestVerC';
[FileNameOutBiVerC] = MakeMeshBidomainVert(Data,N,DomainLimits,FolderName,FileName,'C');


%%
Data = struct;

N = 100;
DomainLimits = [-1 1 -1 1];
FolderName = 'MeshTest/';
FileName = 'QuadriDomain';
[FileNameOutQuadri] = MakeMeshQuadridomain(Data,N,DomainLimits,FolderName,FileName,'P');


N = [20 10];
DomainLimits = [-1 1 -1 1];
FolderName = 'MeshTest/';
FileName = 'QuadriDomainC';
[FileNameOutQuadriC] = MakeMeshQuadridomain(Data,N,DomainLimits,FolderName,FileName,'C');


