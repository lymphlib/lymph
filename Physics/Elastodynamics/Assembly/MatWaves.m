%> @file  MatWaves.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Assembly of the matrices for the elastodinamics problem
%>
%==========================================================================
%> @section classMatWaves Class description
%> @brief  Assembly of the matrices for waves's problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Matrices.Ela     = Matrices for elastic domain
%>                   
%==========================================================================

function [Matrices] = MatWaves(Data, neighbor, femregion)

[Matrices] = MatEla(Data, neighbor, femregion);


