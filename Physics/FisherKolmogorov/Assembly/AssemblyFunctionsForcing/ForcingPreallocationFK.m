%> @file   ForcingPreallocationFK.m
%> @author Mattia Corti
%> @date   16 September 2025
%> @brief Preallocation of the forcing term structures for FK problem.
%>
%==========================================================================
%> @section classForcingPreallocationFK Class description
%==========================================================================
%> @brief           Preallocation of the forcing term structures for FK problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - Vector: Preallocated vector terms
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored 
%> using cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationFK(GenMatrices)
    
    %% Initialization of the forcing terms
    Forcing.Volume.F_loc  = GenMatrices.Vector;

    Forcing.Faces.F_D_loc = GenMatrices.Vector;
    
end
    
