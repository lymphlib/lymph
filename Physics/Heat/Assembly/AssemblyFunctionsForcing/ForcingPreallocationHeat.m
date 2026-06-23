%> @file   ForcingPreallocationHeat.m
%> @author Mattia Corti
%> @date   5 February 2026
%> @brief Preallocation of the vector structure for heat equation.
%>
%==========================================================================
%> @section classForcingPreallocationHeat Class description
%==========================================================================
%> @brief           Preallocation of the vector structure for heat equation 
%> problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - Vector: Preallocated vector terms
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored 
%> using cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationHeat(GenMatrices)

    %% Initialization of the forcing terms
    Forcing.Volume.F_loc  = GenMatrices.Vector;

    Forcing.Faces.F_D_loc = GenMatrices.Vector;
    
end
