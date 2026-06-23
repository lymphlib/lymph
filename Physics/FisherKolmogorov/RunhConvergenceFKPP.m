%> @file  RunhConvergenceFKPP.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Convergence analysis for the Fisher KPP problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergenceFKPP Class description
%==========================================================================
%> @brief          Sequence of run of MainFKPP.m
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
%DataConvTestDir;
DatapAdaptiveWaveFKPP;

%% Initialization error structure
Errors = struct('err_c_L2',[],'err_c_dG',[],'err_Energy',[],'h',[],'NDoF',[]);

for ii = 1:length(Data.meshfileseq)

    Data.name = ['TestFKNeup1_Ref_',num2str(ii)];
    %% Mesh Generation
    if Data.meshfromfile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.foldername,Data.meshfileseq(ii));
    else
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N{ii},Data.domain,Data.foldername,Data.meshfileseq(ii),'P');
    end
    
    %% Main 
    [Error] = MainFKPP(Data,Setup)
    
    %% Errors update
    format shorte

    Errors.err_c_L2     = [Errors.err_c_L2, Error.L2];
    Errors.err_c_dG     = [Errors.err_c_dG, Error.dG];
    Errors.err_Energy   = [Errors.err_Energy, Error.err_Energy];
    Errors.h            = [Errors.h, Error.h];
    Errors.NDoF     = [Errors.NDoF Error.NDoF];


end

%% Plot of the errors

    figure

    cols = lines(3);
    leg = strings(0);
    light = @(c,a) c + a*(1-c);
    c_i= cols(1,:);
    
if Data.adaptivity

    loglog(sqrt(Errors.NDoF),  Errors.err_c_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:));
    grid on
    hold on
    loglog(sqrt(Errors.NDoF), Errors.err_c_dG, '-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:));

    legend("Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
    set(gca,'XScale','log','YScale','log');

    xlabel('$\mathrm{NDoF}^{1/2}$','Interpreter','latex');
    ylabel('$\|u_h - u_{\mathrm{ex}}\|$','Interpreter','latex');

else


    loglog(Errors.h,Errors.h.^Data.degree,'k','LineWidth',1.2)
    hold on
    loglog(Errors.h,Errors.h.^(Data.degree+1),'k--','LineWidth',1.2)
    loglog(Errors.h,Errors.err_c_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:));
    loglog(Errors.h,Errors.err_c_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:));
    loglog(Errors.h,Errors.err_Energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:));

    conv1 = ['$h^', num2str(Data.degree), '$'];
    conv2 = ['$h^', num2str(Data.degree+1), '$'];
    legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
    grid on
    Errors.order_L2 = log(Errors.err_c_L2(1:end-1)./Errors.err_c_L2(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
    Errors.order_dG = log(Errors.err_c_dG(1:end-1)./Errors.err_c_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
    Errors.order_Energy = log(Errors.err_Energy(1:end-1)./Errors.err_Energy(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end))

end