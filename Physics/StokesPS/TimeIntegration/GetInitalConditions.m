%> @file  GetInitalConditions.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Compute initial conditions
%>
%==========================================================================
%> @section classGetInitalConditionsPS Class description
%> @brief  GetInitalConditions
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Matrices    Struct with problem's matrices 
%
%> @retval Uold       Vector with initial conditions
%> @retval Solutions  Struct with initial conditions for saving purposes
%>
%==========================================================================

function [Uold, Solutions] = GetInitalConditions(Data,femregion, Matrices)

[sigma0] = ComputeModalSolutionStokesPS(Data,femregion);

sigma0t = sigma0 * Data.sigma_t_ex{1}(0);
dsigma0t = sigma0 * Data.sigma_dt_ex{1}(0);

sigma0t  = Matrices.Fluid.MPrjP\sigma0t;
dsigma0t = Matrices.Fluid.MPrjP\dsigma0t;


Uold = sigma0t;
Solutions.sigma_h = sigma0t;  
Solutions.dot_sigma_h = dsigma0t;  
