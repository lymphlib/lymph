%> @file  RunpConvergenceFKPP.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Convergence analysis for the Fisher-KPP problem (polynomial refinements)
%>
%==========================================================================
%> @section classDecr Class description
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
DataConvWavesFKPP;

degree_vector = Data.degree;

%% Initialization error structure
Errors = struct('err_c_L2',[],'err_c_dG',[],'err_Energy',[],'p',Data.degree);

%% Parallel pool initialization
if Data.parallel
    delete(gcp('nocreate'));
    parpool(Data.cores);
end

%% Mesh Generation
if Data.meshfromfile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.foldername,Data.meshfileseq);
else
    [Data.meshfile] = MakeMeshMonodomain(Data,Data.N,Data.domain,Data.foldername,Data.meshfileseq,'P');
end

for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

    %% Main
    [Error] = MainFKPP(Data,Setup)

    %% Errors update
    Errors.err_c_L2 = [Errors.err_c_L2, Error.L2];
    Errors.err_c_dG = [Errors.err_c_dG, Error.dG];
    Errors.err_Energy = [Errors.err_Energy, Error.err_Energy];

    format short
    Errors
end

%% Plot of the errors
figure
cols = lines(3);
leg = strings(0);
light = @(c,a) c + a*(1-c);
c_i= cols(1,:);

semilogy(degree_vector,Errors.err_c_L2,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
hold on
semilogy(degree_vector,Errors.err_c_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
semilogy(degree_vector,Errors.err_Energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:))
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on
