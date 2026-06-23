%> @file  RunhConvergenceLaplacian.m
%> @author The Lymph Team
%> @date 5 June 2026
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

%% Initial Simulation Setup
MyPhysicsPath = pwd;
run('../SimulationSetup.m');

%% Input Data - Boundary conditions - Forcing term
%DataConvTestLaph;
DatapAdaptiveConvLap;

Errors.err_L2 = [];
Errors.err_dG = [];
Errors.h = [];
Errors.NDoF = [];

%% Mesh Generation
for ii = 1:4

    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq{ii});
    else
        % Create a new mesh
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N{ii},Data.domain,Data.FolderName,Data.meshfileseq{ii},'P');
    end
    
    
    %% Main     
    [Error] = MainLaplacian(Data,Setup);
    Errors.err_L2   = [Errors.err_L2, Error.L2];
    Errors.err_dG   = [Errors.err_dG, Error.dG];
    Errors.h        = [Errors.h, Error.h];
    Errors.NDoF     = [Errors.NDoF Error.NDoF];

end

%% Plot of the errors
figure
cols = lines(2);
leg = strings(0);
light = @(c,a) c + a*(1-c);
c_i= cols(1,:);

if Data.Adaptivity    

    loglog(sqrt(Errors.NDoF),  Errors.err_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:));
    grid on
    hold on
    loglog(sqrt(Errors.NDoF), Errors.err_dG, '-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:));

    legend("Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
    set(gca,'XScale','log','YScale','log');

    xlabel('$\mathrm{NDoF}^{1/2}$','Interpreter','latex');
    ylabel('$\|u_h - u_{\mathrm{ex}}\|$','Interpreter','latex');

else

    loglog(Errors.h,Errors.h.^Data.degree,'k:','Linewidth',2)
    hold on
    loglog(Errors.h,Errors.h.^(Data.degree+1),'k--','Linewidth',2)
    loglog(Errors.h,Errors.err_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
    loglog(Errors.h,Errors.err_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
    xlabel('h');
    conv1 = ['$h^', num2str(Data.degree), '$'];
    conv2 = ['$h^', num2str(Data.degree+1), '$'];
    legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
    grid on
    Errors.order_L2 = log(Errors.err_L2(1:end-1)./Errors.err_L2(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
    Errors.order_dG = log(Errors.err_dG(1:end-1)./Errors.err_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));

end
