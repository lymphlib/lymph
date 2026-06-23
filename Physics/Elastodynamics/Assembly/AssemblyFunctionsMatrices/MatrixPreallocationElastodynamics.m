%> @file   MatrixPreallocationElastodynamics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   8 May 2026
%> @brief Preallocation of the matrices structure for the elastodynamics problem.
%>
%==========================================================================
%> @section classMatrixPreallocationElastodynamics Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for elastodynamics 
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

function [Matrices] = MatrixPreallocationElastodynamics(GenMatrices)

    %% Initialization of the volume matrices
    Matrices.Volume.V1_loc = GenMatrices.VolMatrix;
    Matrices.Volume.V2_loc = GenMatrices.VolMatrix;
    Matrices.Volume.V3_loc = GenMatrices.VolMatrix;
    Matrices.Volume.V4_loc = GenMatrices.VolMatrix;

    Matrices.Volume.M1_P_rho_loc = GenMatrices.VolMatrix;
    Matrices.Volume.MPrjP_1_loc  = GenMatrices.VolMatrix;
  
    Matrices.Volume.D1_loc = GenMatrices.VolMatrix;
    Matrices.Volume.C1_loc = GenMatrices.VolMatrix;
  
    %% Initialization of the face matrices
    Matrices.Faces.ABC_R1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_R2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_R3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_R4_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.ABC_S1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_S2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_S3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_S4_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.S1_P_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S4_P_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.IT1_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.IT2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.IT3_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.IT4_loc = GenMatrices.FaceMatrix;

end
    