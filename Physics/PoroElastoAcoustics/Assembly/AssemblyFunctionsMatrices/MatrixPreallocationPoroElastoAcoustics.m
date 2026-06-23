%> @file   MatrixPreallocationPoroElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the matrices structure for the poroelastoacoustics problem.
%>
%==========================================================================
%> @section classMatrixPreallocationPoroElastoAcoustics Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for poroelastoacoustics 
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

function [Matrices] = MatrixPreallocationPoroElastoAcoustics(GenMatrices)

    %% Preallocation of the matrices for the elastodynamics (from Physics/Elastodynamics)
    [MatricesEla] = MatrixPreallocationElastodynamics(GenMatrices);
    
    %% Preallocation of the matrices for the poroelasticity
    [MatricesPoro] = MatrixPreallocationPoroelasticity(GenMatrices);

    %% Preallocation of the matrices for the acoustics
    [MatricesAcu] = MatrixPreallocationAcoustics(GenMatrices);

    %% Preallocation of the matrices for the coupling
    [MatricesElaAcu]  = MatrixPreallocationElastoAcoustics(GenMatrices);
    [MatricesPoroAcu] = MatrixPreallocationPoroAcoustics(GenMatrices);
    [MatricesPoroEla] = MatrixPreallocationPoroElastic(GenMatrices);

    %% Creation of final matrices structure
    Matrices.Volume.Ela  = MatricesEla.Volume;
    Matrices.Volume.Poro = MatricesPoro.Volume;
    Matrices.Volume.Acu  = MatricesAcu.Volume;

    Matrices.Faces.Ela  = MatricesEla.Faces;
    Matrices.Faces.Poro = MatricesPoro.Faces;
    Matrices.Faces.Acu  = MatricesAcu.Faces;

    Matrices.Faces.ElaAcu  = MatricesElaAcu.Faces;
    Matrices.Faces.PoroAcu = MatricesPoroAcu.Faces;
    Matrices.Faces.PoroEla = MatricesPoroEla.Faces;


end
    
