%> @file   MatrixPreallocationPoroElastic.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the matrices structure for the poroelasticity-elastodynamics coupling.
%>
%==========================================================================
%> @section classMatrixPreallocationPoroElastic Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for poroelasticity-elastodynamics coupling.
%>
%> @param GenMatrices Struct containing the possible matrices to assemble:
%>                     - VolMatrix: Preallocated matrix for volume terms
%>                     - FaceMatrix: Preallocated matrix for faces terms
%>
%> @retval Matrices  Matrices struct containing local matrices stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Matrices] = MatrixPreallocationPoroElastic(GenMatrices)

    %% Initialization of the coupling matrices
    Matrices.Faces.C1_P_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C2_P_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C3_P_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C4_P_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.C1_E_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C2_E_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C3_E_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.C4_E_loc = GenMatrices.FaceMatrix;

end
    