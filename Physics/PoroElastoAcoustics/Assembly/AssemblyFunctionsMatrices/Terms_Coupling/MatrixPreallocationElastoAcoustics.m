%> @file   MatrixPreallocationElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the matrices structure for the elastodynamics-acoustic coupling.
%>
%==========================================================================
%> @section classMatrixPreallocationElastoAcoustics Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for elastodynamics-acoustic coupling.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - VolMatrix: Preallocated matrix for volume terms
%>                     - FaceMatrix: Preallocated matrix for faces terms
%>
%> @retval Matrices  Matrices struct containing local matrices stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Matrices] = MatrixPreallocationElastoAcoustics(GenMatrices)

    %% Initialization of the coupling matrices
    Matrices.Faces.C1_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C2_loc = GenMatrices.FaceMatrix;

end
    