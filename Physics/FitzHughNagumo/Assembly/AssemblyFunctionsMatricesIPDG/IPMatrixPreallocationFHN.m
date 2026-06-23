%> @file   IPMatrixPreallocationFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief Preallocation of the matrices structures for FHN problem with IP method solver.
%>
%==========================================================================
%> @section classIPMatrixPreallocationFHN Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for FHN 
%> problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - VolMatrix:   Preallocated matrix for volume terms
%>                     - VolMatrix3L: Preallocated matrix for volume terms associated with trilinear terms
%>                     - FaceMatrix:  Preallocated matrix for faces terms
%>
%> @retval Matrices  Matrices struct containing local matrices stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================


function [Matrices] = IPMatrixPreallocationFHN(GenMatrices)
    
    %% Initialization of the volume matrices
    Matrices.Volume.M_prj_loc = GenMatrices.VolMatrix;
    Matrices.Volume.M_u_loc   = GenMatrices.VolMatrix;
    Matrices.Volume.M_w_loc   = GenMatrices.VolMatrix;
    Matrices.Volume.A_loc     = GenMatrices.VolMatrix;

    %% Initialization of the face matrices
    Matrices.Faces.IA_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.SA_loc = GenMatrices.FaceMatrix;

end
    
