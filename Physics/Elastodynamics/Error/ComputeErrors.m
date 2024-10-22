%> @file  ComputeErrors.m
%> @author Ilario Mazzieri, Stefano Bonetti
%> @date 24 July 2024
%> @brief  Compute errors for the elastodynamics problem
%>
%==========================================================================
%> @section classElastodynamicsComputeErrors Class description
%> @brief  Compute errors for waves's problem
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Struct with problem's matrices 
%> @param Solution   Solution vector
%
%> @retval Error     Struct with computed errors
%>
%==========================================================================

function [Error] = ComputeErrors(Data, femregion, Matrices, Solution)


%% Compute modal coefficient of the exact solution
[ue_ex, ve_ex] = ComputeModalSolutionWave(Data,femregion,Data.T);

% Project from nodal to modal representation of the solution coefficients
ue_ex = Matrices.Ela.MPrjP \ ue_ex;
ve_ex = Matrices.Ela.MPrjP \ ve_ex;

%% Definition of error vector
error_ue = ue_ex - Solution(1 : 2*femregion.ndof_e);
error_ve = ve_ex - Solution(2*femregion.ndof_e+1 : end);

%% DG error (displacement)
err_u_dG = error_ue' * Matrices.Ela.DGe * error_ue;    

%% L2 errors
err_u_L2 = error_ue' * Matrices.Ela.MPrjP * error_ue;
err_v_L2 = error_ve' * Matrices.Ela.M_P_rho * error_ve;              

%% Total errors 

err_energy = sqrt(err_v_L2 + err_u_dG); % Energy error
err_u_L2   = sqrt(err_u_L2);
err_u_dG   = sqrt(err_u_dG);
err_v_L2   = sqrt(err_v_L2);

%% Outputs
% Computed errors
Error.err_u_L2   = err_u_L2;
Error.err_u_dG   = err_u_dG;
Error.err_v_L2   = err_v_L2;
Error.err_energy = err_energy;

% Discretization variables
Error.nel = femregion.nel;
Error.h   = Data.h;
Error.p   = Data.degree;
