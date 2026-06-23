%> @file  RunhConvergenceFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 June 2026
%> @brief Convergence analysis for the FHN problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergenceFHN Class description
%==========================================================================
%> @brief          Sequence of run of MainFHN.m
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
DatahConvergenceTestFHN

%% Initialization error structure
Errors = struct('err_c_L2',[],'err_c_dG',[],'h',[]);

for ii = 1:length(Data.meshfileseq)

    Data.name = ['TestFKNeup1_Ref_',num2str(ii)];
    %% Mesh Generation
    if Data.meshfromfile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.foldername,Data.meshfileseq(ii));
    else
        [Data.meshfile] = MakeMeshMonodomain(Data,Data.N(ii),Data.domain,Data.foldername,Data.meshfileseq(ii),'P');
    end
    
    %% Main 
    [Error] = MainFHN(Data,Setup)
    
    %% Errors update
    format shorte

    Errors.err_c_L2     = [Errors.err_c_L2, Error.L2];
    Errors.err_c_dG     = [Errors.err_c_dG, Error.dG];
    Errors.h            = [Errors.h, Error.h];

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
loglog(Errors.h,Errors.err_c_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
loglog(Errors.h,Errors.err_c_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
conv1 = ['$h^', num2str(Data.degree), '$'];
conv2 = ['$h^', num2str(Data.degree+1), '$'];
legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
grid on

%% Computation of convergence rates
Errors.order_L2 = log(Errors.err_c_L2(1:end-1)./Errors.err_c_L2(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end))
Errors.order_dG = log(Errors.err_c_dG(1:end-1)./Errors.err_c_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end))
