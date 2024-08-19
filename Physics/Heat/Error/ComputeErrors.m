%> @file  ComputeErrors.m
%> @author Mattia Corti
%> @date 26 July 2024
%> @brief Compute errors for convergence analysis
%>
%==========================================================================
%> @section classHeatComputeErrors Class description
%==========================================================================
%> @brief Compute errors for convergence analysis
%
%> @param Data        Struct with problem's data
%> @param Matrices    Finite element matrices
%> @param femregion   Structure containing all the information
%> about the finite element approximation
%> @param u_h         Numerical solution of the problem
%> @param time        Time associated to the solution \f$u_h\f$
%>
%> @retval Error      Structure with computed \f$L^2\f$ and dG errors
%>
%==========================================================================

function [Errors] = ComputeErrors(Data, Matrices, femregion, u_h, time)

    %% Compute the exact solutions at the current time
    [u_ex] = EvaluateSolution(Data, femregion, time);

    u_ex = Matrices.M_prj\u_ex;

    %% Compute errors

    error_u = u_h - u_ex;
    err_u_L2 = error_u' * Matrices.M_prj * error_u;
    err_u_dG = error_u' * Matrices.dGA * error_u;

    %% Assign results to the output structs

    Errors.err_u_L2 = full(sqrt(err_u_L2));
    Errors.err_u_dG = full(sqrt(err_u_dG));

    % Discretization variables
    Errors.nel = femregion.nel;
    Errors.h   = Data.h;
    Errors.p   = Data.degree;
    
end
