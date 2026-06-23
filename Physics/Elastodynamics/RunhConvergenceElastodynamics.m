%> @file  RunhConvergenceElastodynamics.m
%> @author The Lymph Team
%> @date 5 June 2026
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


%% Initial Simulation Setup
MyPhysicsPath = pwd;
run('../SimulationSetup.m');

%% Input Data - Boundary conditions - Forcing term
DataConvhTestEla;

Errors.err_u_L2   = [];
Errors.err_u_dG   = [];
Errors.err_v_L2   = [];
Errors.err_energy = [];
Errors.h          = [];

%% Mesh Generation
for ii = 1:4

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq(ii));
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N{ii},Data.domain,Data.FolderName,Data.meshfileseq{ii},'P');
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

cols = lines(3);
leg = strings(0);
light = @(c,a) c + a*(1-c); 
c_i= cols(1,:);

loglog(Errors.h,Errors.h.^Data.degree,'k','LineWidth',1.2)
hold on
loglog(Errors.h,Errors.h.^(Data.degree+1),'k--','LineWidth',1.2)
loglog(Errors.h,Errors.err_u_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:));
loglog(Errors.h,Errors.err_u_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:));
loglog(Errors.h,Errors.err_energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:));

conv1 = ['$h^', num2str(Data.degree), '$'];
conv2 = ['$h^', num2str(Data.degree+1), '$'];
legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
Errors.order_L2 = log(Errors.err_u_L2(1:end-1)./Errors.err_u_L2(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_dG = log(Errors.err_u_dG(1:end-1)./Errors.err_u_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_Energy = log(Errors.err_energy(1:end-1)./Errors.err_energy(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end))

