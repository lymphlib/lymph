%> @file   ForcingPreallocationAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the vector structure for acoustics equation.
%>
%==========================================================================
%> @section classForcingPreallocationAcoustics Class description
%==========================================================================
%> @brief           Preallocation of the vector structure for acoustics equation 
%> problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - Vector: Preallocated vector terms
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored 
%> using cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationAcoustics(GenMatrices)

    %% Initialization of the forcing terms
    Forcing.Volume.F1_loc  = GenMatrices.Vector;

    Forcing.Faces.F1_D_loc = GenMatrices.Vector;
    
end
