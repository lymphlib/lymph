%> @file  RunpConvergenceLaplacian.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Convergence analysis for the Poisson problem (polynomial refinement)
%>
%==========================================================================
%> @section classRunpConvergenceLaplacian Class description
%==========================================================================
%> @brief          Sequence of run of MainLaplacian.m
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
DataConvTestLapp;

Errors.err_L2 = [];
Errors.err_dG = [];
Errors.h = [];

%% Mesh Generation
if Data.MeshFromFile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq);
else
    % Create a new mesh
    [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
end

%p-convergence
for ii = 1:length(Data.degree_vector)

    Data.degree = Data.degree_vector(ii);
    
    %% Main     
    [Error] = MainLaplacian(Data,Setup);
    Errors.err_L2   = [Errors.err_L2, Error.L2];
    Errors.err_dG   = [Errors.err_dG, Error.dG];
    Errors.h        = [Errors.h, Error.h];

end

%% Plot of the errors
figure
cols = lines(2);
leg = strings(0);
light = @(c,a) c + a*(1-c); 
c_i= cols(1,:);

semilogy(Data.degree_vector,Errors.err_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
hold on
semilogy(Data.degree_vector,Errors.err_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
legend("Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
grid on

