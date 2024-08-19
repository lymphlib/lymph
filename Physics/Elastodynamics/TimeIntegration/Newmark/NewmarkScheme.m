%> @file  NewmarkScheme.m
%> @author Ilario Mazzieri, Stefano Bonetti
%> @date 24 July 2024
%> @brief  Newmark time integration scheme
%>
%==========================================================================
%> @section classNewmarkScheme Class description
%> @brief  Newmark time integration scheme
%
%> @param Setup       Struct with problem's setup
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param mesh        Struct with mesh geometry for output
%> @param A,B,C       Matrices and vector for AX''+ BX'+ CX = F
%> @param U_old       Initial vector [X(0), X'(0)]
%
%> @retval U_h        Solution vector (final time)
%>
%==========================================================================

function [U_h] = NewmarkScheme(Setup, Data, femregion, mesh, A, B, C, U_old)

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

clear A B C

% Right-hand side assembly
[F_old, G_old] = ForEla(Data, mesh.neighbor, femregion, Data.t0);

t = 0;
disp(['Starting time: ', num2str(t)]);
disp('------------------------------------------------------')

tic;
counter = 1;

for t = dt : dt : nts*dt
    
    disp(['time: ', num2str(t)]);
    
    [F_new, G_new] = ForEla(Data, mesh.neighbor, femregion, t);
    
    rhs = [dt^2*beta_nm*(F_new + G_new) + dt^2*(1/2-beta_nm)*(F_old + G_old); ...
           gamma_nm*dt*(F_new + G_new)  + dt*(1-gamma_nm)*(F_old + G_old)];

    U_h = A1\(A2*U_old + rhs);
    
    if (mod(counter,Data.VisualizationStep)==0) && (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution)
          PostProcessSolution(Setup, Data, mesh, femregion, counter, U_h,t);
    end
    
    U_old   = U_h;
    F_old   = F_new;
    G_old   = G_new;
    counter = counter + 1;

end