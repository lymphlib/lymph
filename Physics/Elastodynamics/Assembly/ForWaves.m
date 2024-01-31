%> @file  ForWaves.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Assembly of the rhs for the elastodynamics problem
%>
%==========================================================================
%> @section classForWaves Class description
%> @brief  Assembly of the matrices for waves's problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval F       Struct containing the rhs for wave's problem
%(integration wrt the space variable)
%==========================================================================

function [F] = ForWaves(Data, neighbor, femregion)

[F] = ForEla(Data, neighbor, femregion);


