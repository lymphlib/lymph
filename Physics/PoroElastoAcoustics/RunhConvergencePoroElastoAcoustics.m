%> @file  RunhConvergencePoroElastoAcoustics.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Convergence analysis for the coupled poro-elasto-acoustic problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergencePoroElastoAcoustics Class description
%==========================================================================
%> @brief          Sequence of run of RunhConvergencePoroElastoAcoustics.m
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================


%% Initial Simulation Setup
MyPhysicsPath = pwd;
run('../SimulationSetup.m');

% Reuse of the assembly routines of elastodynamics (excluded faces assembly) to avoid code duplication
addpath(genpath('..\Elastodynamics\Assembly\'));

%% Input Data - Boundary conditions - Forcing term
% DataConvTestAcu;
% DataConvTestEla;
% DataConvTestPoro;
 DataConvTestPoroAcu;
% DataConvTestElaAcu
% DataConvTestPoroEla

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
        if (Data.NPhys == 1)
            Data.meshfile = MakeMeshMonodomain(Data,Data.N(ii),Data.domain,Data.FolderName,Data.meshfileseq,'P');
        elseif (Data.NPhys == 2)
            Data.meshfile = MakeMeshBidomainVert(Data,Data.N(ii),Data.domain,Data.FolderName,Data.meshfileseq,'P');
        end
    end

    %% Main
    [Error] = MainPoroAcuEla(Data,Setup);
    Errors.err_L2_v   = [Errors.err_L2_v,   Error.error_L2_v];
    Errors.err_L2_d   = [Errors.err_L2_d,   Error.error_L2_d];
    Errors.err_dG     = [Errors.err_dG,     Error.error_dG];
    Errors.err_Energy = [Errors.err_Energy, Error.error_Energy];
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
loglog(Errors.h,Errors.err_L2_d,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:));
loglog(Errors.h,Errors.err_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:));
loglog(Errors.h,Errors.err_Energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:));

conv1 = ['$h^', num2str(Data.degree), '$'];
conv2 = ['$h^', num2str(Data.degree+1), '$'];
legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
Errors.order_L2 = log(Errors.err_L2_d(1:end-1)./Errors.err_L2_d(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_dG = log(Errors.err_dG(1:end-1)./Errors.err_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_Energy = log(Errors.err_Energy(1:end-1)./Errors.err_Energy(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end))
