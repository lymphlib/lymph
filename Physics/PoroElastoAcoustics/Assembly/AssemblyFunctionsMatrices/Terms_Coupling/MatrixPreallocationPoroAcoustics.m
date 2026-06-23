%> @file   MatrixPreallocationPoroAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   8 May 2026
%> @brief Preallocation of the matrices structure for the poroelasticity-acoustics coupling.
%>
%==========================================================================
%> @section classMatrixPreallocationPoroAcoustics Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for poroelasticity-acoustics coupling.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - VolMatrix: Preallocated matrix for volume terms
%>                     - FaceMatrix: Preallocated matrix for faces terms
%>
%> @retval Matrices  Matrices struct containing local matrices stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Matrices] = MatrixPreallocationPoroAcoustics(GenMatrices)

    %% Initialization of the coupling matrices
    Matrices.Faces.C1_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C2_loc = GenMatrices.FaceMatrix;

end
    