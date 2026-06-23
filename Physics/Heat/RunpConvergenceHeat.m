%> @file  RunpConvergenceHeat.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Convergence analysis for the heat equation (polynomial refinements)
%>
%==========================================================================
%> @section classRunpConvergenceHeat Class description
%==========================================================================
%> @brief          Sequence of run of MainHeat.m
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
DatapConvergenceTest;

degree_vector = Data.degree;

Errors.err_u_L2 = [];
Errors.err_u_dG = [];
Errors.err_Energy = [];
Errors.h = [];

%% Mesh Generation
if Data.MeshFromFile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq);
else
    [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
end
    
%% p convergence
for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);    
    
    %% Main 
    [Error] = MainHeat(Data,Setup);
    Errors.err_u_L2 = [Errors.err_u_L2, Error.L2];
    Errors.err_u_dG = [Errors.err_u_dG, Error.dG];
    Errors.err_Energy = [Errors.err_Energy, Error.err_Energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
cols = lines(3);
leg = strings(0);
light = @(c,a) c + a*(1-c); 
c_i= cols(1,:);

semilogy(degree_vector,Errors.err_u_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
hold on
semilogy(degree_vector,Errors.err_u_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
semilogy(degree_vector,Errors.err_Energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:))
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
