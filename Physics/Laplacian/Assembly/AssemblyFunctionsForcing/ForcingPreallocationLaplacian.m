%> @file   ForcingPreallocationLaplacian.m
%> @author Mattia Corti
%> @date   11 June 2025
%> @brief Preallocation of the vector structure for Laplacian problem.
%>
%==========================================================================
%> @section classForcingPreallocationLaplacian Class description
%==========================================================================
%> @brief           Preallocation of the vector structure for Laplacian 
%> problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - Vector: Preallocated vector terms
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored 
%> using cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationLaplacian(GenMatrices)

    %% Initialization of the forcing terms
    Forcing.Volume.F_loc  = GenMatrices.Vector;

    Forcing.Faces.F_D_loc = GenMatrices.Vector;
    
end
