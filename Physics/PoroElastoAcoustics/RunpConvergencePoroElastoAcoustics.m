%> @file  RunpConvergencePoroElastoAcoustics.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Convergence analysis for the coupled poro-elasto-acoustic problem (p refinements)
%>
%==========================================================================
%> @section classRunpConvergencePoroElastoAcoustics Class description
%==========================================================================
%> @brief          Sequence of run of RunpConvergencePoroElastoAcoustics.m
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

degree_vector = [1 2 3 4 5];

Errors.err_L2_d = [];
Errors.err_L2_v = [];
Errors.err_dG = [];
Errors.err_Energy = [];
Errors.h = [];


%% Mesh Generation

if Data.MeshFromFile
    % Load an existing mesh
    Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq(1));
else
    % Create a new mesh
    if (Data.NPhys == 1)
        Data.meshfile = MakeMeshMonodomain(Data,Data.N(1),Data.domain,Data.FolderName,Data.meshfileseq,'P');
    elseif (Data.NPhys == 2)
        Data.meshfile = MakeMeshBidomainVert(Data,Data.N(1),Data.domain,Data.FolderName,Data.meshfileseq,'P');
    end
end

for ii = 1:length(degree_vector)

    Data.degree = degree_vector(ii);

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

semilogy(degree_vector,Errors.err_L2_d,'-o', 'LineWidth',1.4,'Color',cols(1,:),'MarkerFaceColor',cols(1,:))
hold on
semilogy(degree_vector,Errors.err_dG,'-o', 'LineWidth',1.4,'Color',cols(2,:),'MarkerFaceColor',cols(2,:))
semilogy(degree_vector,Errors.err_Energy,'-o', 'LineWidth',1.4,'Color',cols(3,:),'MarkerFaceColor',cols(3,:))
legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
grid on


% figure
% semilogy(degree_vector,Errors.err_L2_d,'g')
% hold on
% semilogy(degree_vector,Errors.err_dG,'r')
% semilogy(degree_vector,Errors.err_Energy,'b')
% legend("Error $L^2$-norm", "Error DG-norm", "Error energy-norm","Interpreter","latex")
% grid on
