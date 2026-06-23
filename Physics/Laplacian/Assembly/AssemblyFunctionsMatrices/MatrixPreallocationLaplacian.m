%> @file   MatrixPreallocationLaplacian.m
%> @author Mattia Corti
%> @date   11 June 2025
%> @brief Preallocation of the matrices structure for the Laplacian problem.
%>
%==========================================================================
%> @section classMatrixPreallocationLaplacian Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for Laplacian 
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

function [Matrices] = MatrixPreallocationLaplacian(GenMatrices)

    %% Initialization of the volume matrices
    Matrices.Volume.Mprj_loc = GenMatrices.VolMatrix;
    Matrices.Volume.A_loc    = GenMatrices.VolMatrix;

    %% Initialization of the face matrices
    Matrices.Faces.IA_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.SA_loc = GenMatrices.FaceMatrix;
   
end
    