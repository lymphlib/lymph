%> @file   ForcingPreallocationElastodynamics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   8 May 2026
%> @brief Preallocation of the vector structure for Elastodynamics equation.
%>
%==========================================================================
%> @section classForcingPreallocationElastodynamics Class description
%==========================================================================
%> @brief           Preallocation of the vector structure for Elastodynamics equation 
%> problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - Vector: Preallocated vector terms
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored 
%> using cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationElastodynamics(GenMatrices)

    %% Initialization of the forcing terms
    Forcing.Volume.F1_loc  = GenMatrices.Vector;
    Forcing.Volume.F2_loc  = GenMatrices.Vector;
    Forcing.Volume.G1_loc  = GenMatrices.Vector;
    Forcing.Volume.G2_loc  = GenMatrices.Vector;

    Forcing.Faces.F1_D_loc = GenMatrices.Vector;
    Forcing.Faces.F2_D_loc = GenMatrices.Vector;
    
end
