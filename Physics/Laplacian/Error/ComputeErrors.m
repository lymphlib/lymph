%> @file  ComputeErrors.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 26 July 2024
%> @brief Compute errors for convergence analysis
%>
%==========================================================================
%> @section classLaplacianComputeErrors Class description
%==========================================================================
%> @brief Compute errors for convergence analysis
%
%> @param Data  Struct with problem's data
%> @param femregion   Structure containing all the information 
%>                    about the finite element approximation
%> @param Matrices    Finite element matrices 
%> @param U           Problem's solution
%
%> @retval Errror     Structure with computed L2 and dG errors 
%>
%==========================================================================
function [Error] = ComputeErrors(Data, femregion, Matrices, U)

%% Compute modal coefficient of the exact solution
[uex_modal] = ComputeModalSolution(Data, femregion);
uex = Matrices.Mprj\uex_modal;

%% Definition of error vector
error_u = uex - U;

%% DG error
error_dGu = sqrt(error_u' * Matrices.dGA * error_u);

%% L2 error
error_L2u = sqrt(error_u' * (Matrices.Mprj) * error_u);

%% Outputs

Error.nel = femregion.nel;
Error.h   = Data.h;
Error.p   = Data.degree;

Error.L2 = error_L2u;
Error.dG = error_dGu;

