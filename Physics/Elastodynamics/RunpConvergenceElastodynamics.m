%> @file  RunpConvergenceElastodynamics.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Convergence analysis for the Elastodynamics problem (polynomial refinements)
%>
%==========================================================================
%> @section classRunpConvergenceElastodynamics Class description
%==========================================================================
%> @brief          Sequence of run of MainEla.m
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
DataConvpTestEla;

degree_vector = [1 2 3 4 5];

Errors.err_u_L2   = [];
Errors.err_u_dG   = [];
Errors.err_v_L2   = [];
Errors.err_energy = [];
Errors.h          = [];

%% Mesh Generation
for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq);
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.FolderName,Data.meshfileseq,'P');
    end
    
    
    %% Main     

    [Error] = MainEla(Data,Setup);
    Errors.err_u_L2   = [Errors.err_u_L2,   Error.err_u_L2];
    Errors.err_u_dG   = [Errors.err_u_dG,   Error.err_u_dG];
    Errors.err_v_L2   = [Errors.err_v_L2,   Error.err_v_L2];
    Errors.err_energy = [Errors.err_energy, Error.err_energy];
    Errors.h = [Errors.h, Error.h];

end

%% Plot of the errors
figure
semilogy(degree_vector,Errors.err_u_L2,'g')
hold on
semilogy(degree_vector,Errors.err_u_dG,'r')
semilogy(degree_vector,Errors.err_energy,'b')
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on

%% Plot of the errors
figure
cols = lines(3);
leg = strings(0);
light = @(c,a) c + a*(1-c); 
c_i= cols(1,:);

semilogy(degree_vector,Errors.err_u_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
hold on
semilogy(degree_vector,Errors.err_u_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
semilogy(degree_vector,Errors.err_energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:))
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
