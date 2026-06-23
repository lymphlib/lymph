%> @file  RunNDoFConvergencep.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Convergence analysis for the FKPP problem (mesh refinements)
%>
%==========================================================================
%> @section classRunNDoFConvergencepFKPP Class description
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
DataComparep

Errors.Adaptive = {};
Errors.Uniform = {};

for pp = 1: length(Data.l_val)

    Errors.Adaptive{pp}.err_L2 = [];
    Errors.Adaptive{pp}.err_dG = [];
    Errors.Adaptive{pp}.h = [];
    Errors.Adaptive{pp}.NDoF = [];

    Errors.Uniform{pp}.err_L2 = [];
    Errors.Uniform{pp}.err_dG = [];
    Errors.Uniform{pp}.h = [];
    Errors.Uniform{pp}.NDoF = [];

    % Construciton of the Adaptivity function
    Data.degree  = pp;
    Data.maxDegree  = pp;
    Data.AdaptFunc  = @(tau_r) floor(1 + (Data.maxDegree-1)*tanh(5*tau_r)); 

    %% Mesh Generation
    for ii = 1:length(Data.meshfileseq)


        if Data.MeshFromFile
            % Load an existing mesh
            Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq{ii});
        else
            % Create a new mesh
            [Data.meshfile] = MakeMeshMonodomain(Data,Data.N{ii},Data.domain,Data.FolderName,Data.meshfileseq{ii},'P');
        end

        %% Main
        Data.Adaptivity = true;
        [Error] = MainFKPP(Data,Setup);

        Errors.Adaptive{pp}.err_L2   = [Errors.Adaptive{pp}.err_L2, Error.L2];
        Errors.Adaptive{pp}.err_dG   = [Errors.Adaptive{pp}.err_dG, Error.dG];
        Errors.Adaptive{pp}.h        = [Errors.Adaptive{pp}.h, Error.h];
        Errors.Adaptive{pp}.NDoF     = [Errors.Adaptive{pp}.NDoF Error.NDoF];

        Data.Adaptivity = false;
        [Error] = MainFKPP(Data,Setup);

        Errors.Uniform{pp}.err_L2   = [Errors.Uniform{pp}.err_L2, Error.L2];
        Errors.Uniform{pp}.err_dG   = [Errors.Uniform{pp}.err_dG, Error.dG];
        Errors.Uniform{pp}.h        = [Errors.Uniform{pp}.h, Error.h];
        Errors.Uniform{pp}.NDoF     = [Errors.Uniform{pp}.NDoF Error.NDoF];

    end

end


%% Plot NDoF
figure

cols = lines(length(Data.l_val));
leg = strings(0);
light = @(c,a) c + a*(1-c);   % light colors

subplot(1,2,1)
for ii=1:length(Data.l_val)
    c_i= cols(ii,:);
    c_iL = light(c_i,0.55);
    loglog(sqrt(Errors.Adaptive{ii}.NDoF),  Errors.Adaptive{ii}.err_L2,'-o', 'LineWidth',1.4,'Color',c_i,'MarkerFaceColor',c_i);
    leg(end+1) = "$p=" + Data.l_val(ii) + "${ ad.}";
    hold on
    loglog(sqrt(Errors.Uniform{ii}.NDoF),  Errors.Uniform{ii}.err_L2,'-o', 'LineWidth',1.4,'Color',c_iL,'MarkerFaceColor',c_iL);
    leg(end+1) = "$p=" + Data.l_val(ii) + "${ un.}";
end
grid on
legend(leg,"Interpreter","latex")
set(gca,'XScale','log','YScale','log');
xlabel('$\mathrm{NDoF}^{1/2}$','Interpreter','latex');
title('$\|u_h - u_{\mathrm{ex}}\|_{L^2}$','Interpreter','latex');

subplot(1,2,2)
for ii=1:length(Data.l_val)
    c_i= cols(ii,:);
    c_iL = light(c_i,0.55);
    loglog(sqrt(Errors.Adaptive{ii}.NDoF),  Errors.Adaptive{ii}.err_dG,'-o', 'LineWidth',1.4,'Color',c_i,'MarkerFaceColor',c_i);
    leg(end+1) = "$p=" + Data.l_val(ii) + "${ ad.}";
    hold on
    loglog(sqrt(Errors.Uniform{ii}.NDoF),  Errors.Uniform{ii}.err_dG,'-o', 'LineWidth',1.4,'Color',c_iL,'MarkerFaceColor',c_iL);
    leg(end+1) = "$p=" + Data.l_val(ii) + "${ un.}";
end
grid on
legend(leg,"Interpreter","latex")
set(gca,'XScale','log','YScale','log');
xlabel('$\mathrm{NDoF}^{1/2}$','Interpreter','latex');
title('$\|u_h - u_{\mathrm{ex}}\|_{dG}$','Interpreter','latex');

