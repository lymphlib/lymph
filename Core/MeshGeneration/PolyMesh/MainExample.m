%% simple example of mesh generation

clear; clc; close all;
addpath Polymesher

%% 
N = 100;
DomainLimits = [-1 4 -1 1];
Circle.Radius = 0.2;
Circle.Center = [0 0];
FolderName = 'MeshTest/';
FileName = 'FlowAroundCilinder';
[FileNameOutCircle] = MakeMeshFlowAroundCilinder(N,DomainLimits,Circle,FolderName,FileName);



%%
N = 100;
DomainLimits = [0 1 0 1];
FolderName = 'MeshTest/';
FileName = 'MonoDomainTest';
[FileNameOutMono] = MakeMeshMonodomain(N,DomainLimits,FolderName,FileName,'P');

N = [10 10];
DomainLimits = [0 1 0 1];
FolderName = 'MeshTest/';
FileName = 'Quad01_quad_ref0';
[FileNameOutMonoC] = MakeMeshMonodomain(N,DomainLimits,FolderName,FileName,'C');

%%
N = 400;
DomainLimits = [-2400 2400 0 4800];
FolderName = 'MeshTest/';
FileName = 'BiDomainTestVer';
[FileNameOutBiVer] = MakeMeshBidomainVert(N,DomainLimits,FolderName,FileName,'P');

N = [20 20];
DomainLimits = [-2400 2400 0 4800];
FolderName = 'MeshTest/';
FileName = 'BiDomainTestVerC';
[FileNameOutBiVerC] = MakeMeshBidomainVert(N,DomainLimits,FolderName,FileName,'C');


%%
N = [4 4];
DomainLimits = [0 4800 -2400 2400];
FolderName = 'MeshTest/';
FileName = 'BiDomainTestHorC';
[FileNameOutBiHorC] = MakeMeshBidomainHor(N,DomainLimits,FolderName,FileName,'C');

N = 3600;
DomainLimits = [0 4800 -2400 2400];
FolderName = 'MeshTest/';
FileName = 'BiDomainTestHor';
[FileNameOutBiHorP] = MakeMeshBidomainHor(N,DomainLimits,FolderName,FileName,'P');



%%
N = 100;
DomainLimits = [-1 1 -1 1];
FolderName = 'MeshTest/';
FileName = 'QuadriDomain';
[FileNameOutQuadri] = MakeMeshQuadridomain(N,DomainLimits,FolderName,FileName,'P');


N = [20 10];
DomainLimits = [-1 1 -1 1];
FolderName = 'MeshTest/';
FileName = 'QuadriDomainC';
[FileNameOutQuadriC] = MakeMeshQuadridomain(N,DomainLimits,FolderName,FileName,'C');


