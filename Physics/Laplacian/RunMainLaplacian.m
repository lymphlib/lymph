%> @file  RunMainLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
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

%% Import lymph and paths of folders related to this problem
run("../ImportLymphPaths.m")
MyPhysicsPath = pwd;
addpath(genpath(fullfile(MyPhysicsPath,'Assembly')));
addpath(genpath(fullfile(MyPhysicsPath,'InputData')));
addpath(genpath(fullfile(MyPhysicsPath,'MainFunctions')));
addpath(genpath(fullfile(MyPhysicsPath,'Error')));
addpath(genpath(fullfile(MyPhysicsPath,'PostProcessing')));

%% Simulation - Setup
run("../RunSetup.m")

%% Input Data - Boundary conditions - Forcing term
domainType = 'Monodomain'		% choose between 'Monodomain' and 'DoubleCircle'
if strcmp(domainType, 'Monodomain')
	DataTestLap;
elseif strcmp(domainType, 'DoubleCircle')
	DataTestDCLap;
else
	error(['Unknown domainType = ', domainType, '\n\tChoose between ''Monodomain'' and ''DoubleCircle'''])
end


%% Mesh Generation

if Data.MeshFromFile
    % Load existing mesh
    Data.meshfile = fullfile(Data.FolderName, Data.meshfileseq);
else
    % Create a new mesh
	if strcmp(domainType, 'Monodomain')
	    Data.meshfile = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P','laplacian');
	elseif strcmp(domainType, 'DoubleCircle')
		Data.meshfile = MakeMeshDoubleCircle(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'laplacian');
	else
		error(['Unknown domainType = ', domainType, '\n\tChoose between ''Monodomain'' and ''DoubleCircle'''])
	end
end



%% Main
[Error] = MainLaplacian(Data,Setup);
