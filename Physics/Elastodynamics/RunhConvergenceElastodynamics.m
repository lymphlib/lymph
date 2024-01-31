%> @file  RunhConvergenceElastodynamics.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Convergence analysis for the Elastodynamics problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergenceElastodynamics Class description
%==========================================================================
%> @brief          Sequence of run of MainEla.m
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
DataConvTestEla;

Errors.err_L2_d = [];
Errors.err_L2_v = [];
Errors.err_dG = [];
Errors.err_Energy = [];
Errors.h = [];

%% Mesh Generation
for ii = 1:4

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq(ii));
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq{ii},'P','laplacian');
    end
    
    
    %% Main     
    [Error] = MainEla(Data,Setup);
    Errors.err_L2_v   = [Errors.err_L2_v,   Error.error_L2_v];
    Errors.err_L2_d   = [Errors.err_L2_d,   Error.error_L2_d];
    Errors.err_dG     = [Errors.err_dG,     Error.error_dG];
    Errors.err_Energy = [Errors.err_Energy, Error.error_Energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
loglog(Errors.h,Errors.h.^Data.degree,'k')
hold on
loglog(Errors.h,Errors.h.^(Data.degree+1),'k--')
loglog(Errors.h,Errors.err_L2_d,'g')
loglog(Errors.h,Errors.err_dG,'r')
loglog(Errors.h,Errors.err_Energy,'b')
conv1 = ['$h^', num2str(Data.degree), '$'];
conv2 = ['$h^', num2str(Data.degree+1), '$'];
legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
Errors.order_L2 = log(Errors.err_L2_d(1:end-1)./Errors.err_L2_d(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_dG = log(Errors.err_dG(1:end-1)./Errors.err_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_Energy = log(Errors.err_Energy(1:end-1)./Errors.err_Energy(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
