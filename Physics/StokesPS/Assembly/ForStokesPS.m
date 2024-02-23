%> @file  ForStokesPS.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Assembly of the rhs for the Stokes problem in pseudo-stress
%formulation
%>
%==========================================================================
%> @section classForStokesPS Class description
%> @brief  Assembly of the rhs for the Stokes problem in pseudo-stress
%formulation
%>
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%>
%> @retval F       Struct containing the rhs for the Stokes problem
%(integration wrt the space variable)
%==========================================================================

function [F] = ForStokesPS(Data, neighbor, femregion)

fprintf('  Fluid media \n');

[F] = ForFluid(Data, neighbor, femregion);



