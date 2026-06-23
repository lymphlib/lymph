%> @file   ForcingPreallocationFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief Preallocation of the forcing term structures for FHN problem.
%>
%==========================================================================
%> @section classForcingPreallocationFHN Class description
%==========================================================================
%> @brief           Preallocation of the forcing term structures for FHN problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - Vector: Preallocated vector terms
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored 
%> using cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationFHN(GenMatrices)
    
    %% Initialization of the forcing terms
    Forcing.Volume.F_loc     = GenMatrices.Vector;
    Forcing.Volume.Iext_loc  = GenMatrices.Vector;
    Forcing.Volume.G_loc     = GenMatrices.Vector;
    
    Forcing.Faces.F_D_loc = GenMatrices.Vector;
    
end
    
