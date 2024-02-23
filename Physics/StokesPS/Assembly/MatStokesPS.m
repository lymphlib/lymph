%> @file  MatStokesPS.m
%> @author Ilario Mazzieri
%> @date 13 February 2024
%> @brief  Assembly of the matrices for the Stokes problem in pseudo-stress
%formulation
%>
%==========================================================================
%> @section classMatStokesPS Class description
%> @brief  Assembly of the matrices for the Stokes problem
%>
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%>
%> @retval Matrices  Matrices.Poro      = Matrices for poroelastic domain
%>                   Matrices.Fluid     = Matrices for fluid domain
%>                   Matrices.PoroFluid = Matrices for poro-fluid coupling
%>
%==========================================================================

function [Matrices] = MatStokesPS(Data, neighbor, femregion)

fprintf('  Fluid media \n');
[Matrices] = MatFluid(Data, neighbor, femregion);


