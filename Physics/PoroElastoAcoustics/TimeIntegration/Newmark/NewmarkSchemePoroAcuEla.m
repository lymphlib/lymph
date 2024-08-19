%> @file  NewmarkSchemePoroAcuEla.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 13 August 2024
%> @brief  Newmark time integration scheme
%>
%==========================================================================
%> @section classNewmarkSchemePoroAcuEla Class description
%> @brief  Newmark time integration scheme
%
%> @param Setup       Struct with problem's setup
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param mesh        Struct with mesh geometry for output
%> @param A,B,C,F     Matrices and vector for AX''+ BX'+ CX = F
%> @param Uold        Initial vector [X(0), X'(0)]
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function [Solutions] = NewmarkSchemePoroAcuEla(Setup, Data, femregion, mesh, A, B, C, F, Uold)

% Time discretization parameters
dt  = Data.dt;
nts = floor(Data.T/dt);

% Newmark parameter
beta_nm  = Data.BetaNM;
gamma_nm = Data.GammaNM;

A1 = [A + dt^2*beta_nm*C, dt^2*beta_nm*B; ...
    gamma_nm*dt*C,      A + gamma_nm*dt*B];

A2 = [A-dt^2*(1/2-beta_nm)*C, dt*A - dt^2*(1/2-beta_nm)*B; ...
    -dt*(1-gamma_nm)*C,     A - dt*(1-gamma_nm)*B];

%clear A B C

Fold = [F.f_p * Data.source_up_t{1}(0)  + F.j_p * Data.source_upd_t{1}(0) + F.f_p_diri * Data.DirBCPoro_up_t{1}(0); ...
    F.g_p * Data.source_wp_t{1}(0)  + F.h_p * Data.source_wpd_t{1}(0) + F.g_p_diri * Data.DirBCPoro_wp_t{1}(0); ...
    F.f_a * Data.source_phi_t{1}(0) + F.f_a_diri * Data.DirBCAcu_t{1}(0); ...
    F.f_e * Data.source_ue_t{1}(0)  + F.g_e * Data.source_ued_t{1}(0) + F.f_e_diri * Data.DirBCEla_t{1}(0)];


t = 0;
disp(['Starting time: ', num2str(t)]);
disp('------------------------------------------------------')

tic;
counter = 1;

for t = dt : dt : nts*dt

    disp(['time: ', num2str(t)]);

    Fnew = [F.f_p * Data.source_up_t{1}(t) + F.j_p * Data.source_upd_t{1}(t) + F.f_p_diri * Data.DirBCPoro_up_t{1}(t); ...
        F.g_p * Data.source_wp_t{1}(t) + F.h_p * Data.source_wpd_t{1}(t) + + F.g_p_diri * Data.DirBCPoro_wp_t{1}(t); ...
        F.f_a * Data.source_phi_t{1}(t) + F.f_a_diri * Data.DirBCAcu_t{1}(t); ...
        F.f_e * Data.source_ue_t{1}(t) + F.g_e * Data.source_ued_t{1}(t) + F.f_e_diri * Data.DirBCEla_t{1}(t)];

    rhs = [dt^2*beta_nm*Fnew + dt^2*(1/2-beta_nm)*Fold; ...
        gamma_nm*dt*Fnew  + dt*(1-gamma_nm)*Fold];

    Uh = A1\(A2*Uold + rhs);

    if (mod(counter,Data.VisualizationStep)==0) && (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution)
        % Construct solution struct
        [Solutions] = SaveSolutionPoroAcuEla(Uh,femregion);

        % Postprocess solution
        PostProcessSolution(Setup, Data, mesh, femregion, counter, Solutions,t);
    end

    Uold    = Uh;
    Fold    = Fnew;
    counter = counter + 1;

end

[Solutions] = SaveSolutionPoroAcuEla(Uh,femregion);

