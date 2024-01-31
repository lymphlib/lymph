%> @file  ComputeErrors.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Compute errors for the elastodynamics problem
%>
%==========================================================================
%> @section classElastodynamicsComputeErrors Class description
%> @brief  Compute errors for waves's problem
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Struct with problem's matrices 
%> @param Solutions  Struct with solution vectors 
%
%> @retval Error     Struct with computed errors
%>
%==========================================================================

function [Error] = ComputeErrors(Data, femregion, Matrices, Solutions)


%% Compute modal coefficient of the exact solution
[ue0] = ComputeModalSolutionWave(Data,femregion);

T = Data.T;
% Solution
ue_ex   = ue0   * Data.ue_t_ex{1}(T);
% Derivative
dot_ue_ex   = ue0   * Data.due_t_ex{1}(T);

%% Switch from nodal to modal representation of the solution coefficients
ue_ex     = Matrices.Ela.MPrjP\ue_ex;
dot_ue_ex = Matrices.Ela.MPrjP\dot_ue_ex;

%% Definition of error vector
error_ue       = ue_ex  - Solutions.ue_h;
error_dot_ue   = dot_ue_ex  - Solutions.dot_ue_h;

%% DG errors
error_dG = error_ue'    * Matrices.Ela.DGe * error_ue;    

%% L2 velocity errors
error_L2_vel   = error_dot_ue'   * (    Matrices.Ela.M_P_rho      ) * error_dot_ue;

%% L2 displacement errors                
error_L2   = error_ue'   * Matrices.Ela.MPrjP   * error_ue;

%% Total errors 

% Energy error
error_Energy = sqrt(error_L2_vel + error_dG);
error_L2_v   = sqrt(error_L2_vel);
error_L2_d   = sqrt(error_L2);


%% Outputs
Error.nel = femregion.nel;
Error.h   = Data.h;
Error.p   = Data.degree;

Error.error_L2_v        = error_L2_v;
Error.error_L2_d        = error_L2_d;
Error.error_Energy      = error_Energy;
Error.error_dG          = sqrt(error_dG);
Error.error_L2          = error_L2;
Error.error_L2_vel      = error_L2_vel;


