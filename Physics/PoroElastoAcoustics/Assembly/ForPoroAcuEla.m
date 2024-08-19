%> @file  ForPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 26 June 2024
%> @brief  Assembly of the rhs for the poro-acoustic-elastic problem
%>
%==========================================================================
%> @section classForPoroAcuEla Class description
%> @brief  Assembly of the matrices for waves's problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval F       Struct containing the rhs for wave's problem
%(integration wrt the space variable)
%==========================================================================

function [F] = ForPoroAcuEla(Data, neighbor, femregion)
%% RIGHT HAND SIDE
% Porous media
disp('Rhs for Porous media ... ');
[F] = ForPoro(Data, neighbor, femregion);
disp('Done')
disp('------------------------------------------------------')

% Acoustic media
disp('Rhs for Acoustic media ... ');
[F] = ForAcu(Data, neighbor, femregion, F);
disp('Done')
disp('------------------------------------------------------')

% Elastic media
disp('Rhs for Elastic media ... ');
[F] = ForEla(Data, neighbor, femregion, F);
disp('Done')
disp('------------------------------------------------------')


