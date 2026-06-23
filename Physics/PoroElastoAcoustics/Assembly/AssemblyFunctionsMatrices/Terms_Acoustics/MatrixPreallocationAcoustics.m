%> @file   MatrixPreallocationAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the matrices structure for the acoustics problem.
%>
%==========================================================================
%> @section classMatrixPreallocationAcoustics Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for acoustics 
%> problem.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - VolMatrix: Preallocated matrix for volume terms
%>                     - FaceMatrix: Preallocated matrix for faces terms
%>
%> @retval Matrices  Matrices struct containing local matrices stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Matrices] = MatrixPreallocationAcoustics(GenMatrices)

    %% Initialization of the volume matrices
    Matrices.Volume.M_P_A_loc = GenMatrices.VolMatrix;
    Matrices.Volume.MPrjA_loc = GenMatrices.VolMatrix;
    Matrices.Volume.W_loc     = GenMatrices.VolMatrix;

    %% Initialization of the face matrices
    Matrices.Faces.S_A_loc  = GenMatrices.FaceMatrix;
    Matrices.Faces.IT_A_loc = GenMatrices.FaceMatrix;

end
    