%> @file  GetInitalConditions.m
%> @author Mattia Corti
%> @date 26 July 2024
%> @brief Compute initial condition
%>
%==========================================================================
%> @section classHeatGetInitalConditions Class description
%> @brief  GetInitalConditions
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Matrices    Struct with problem's matrices 
%
%> @retval Uold       Vector with initial conditions
%>
%==========================================================================

function [u_old] = GetInitalConditions(Data, femregion, Matrices)

    [u_old]  = EvaluateSolution(Data, femregion, Data.t0);

    % Projection for modal coordinates
    u_old  = Matrices.M_prj\u_old;

end
