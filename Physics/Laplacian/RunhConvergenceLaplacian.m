%> @file  RunhConvergenceLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Convergence analysis for the Poisson problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergenceLaplacian Class description
%==========================================================================
%> @brief          Sequence of run of MainLaplacian.m
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
DataConvTestLap;

Errors.err_L2 = [];
Errors.err_dG = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:4

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq{ii});
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N{ii},Data.domain,Data.FolderName,Data.meshfileseq{ii},'P','laplacian');
    end
    
    
    %% Main     
    [Error] = MainLaplacian(Data,Setup);
    Errors.err_L2   = [Errors.err_L2, Error.L2];
    Errors.err_dG   = [Errors.err_dG, Error.dG];
    Errors.h        = [Errors.h, Error.h];

end

%% Plot of the errors
figure
loglog(Errors.h,Errors.h.^Data.degree,'k:','Linewidth',2)
hold on
loglog(Errors.h,Errors.h.^(Data.degree+1),'k--','Linewidth',2)
loglog(Errors.h,Errors.err_L2,'g','Linewidth',2)
loglog(Errors.h,Errors.err_dG,'r','Linewidth',2)
xlabel('h');
conv1 = ['$h^', num2str(Data.degree), '$'];
conv2 = ['$h^', num2str(Data.degree+1), '$'];
legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
grid on
Errors.order_L2 = log(Errors.err_L2(1:end-1)./Errors.err_L2(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_dG = log(Errors.err_dG(1:end-1)./Errors.err_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
