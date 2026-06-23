%> @file   IPMatrixPreallocationFK.m
%> @author Mattia Corti
%> @date   16 September 2025
%> @brief Preallocation of the matrices structures for Fisher-KPP problem with IP method solver.
%>
%==========================================================================
%> @section classIPMatrixPreallocationFK Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for Fisher-KPP 
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


function [Matrices] = IPMatrixPreallocationFK(GenMatrices)
    
    %% Initialization of the volume matrices
    Matrices.Volume.M_prj_loc = GenMatrices.VolMatrix;
    Matrices.Volume.M_loc     = GenMatrices.VolMatrix;
    Matrices.Volume.A_loc     = GenMatrices.VolMatrix;

    %% Initialization of the face matrices
    Matrices.Faces.IA_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.SA_loc = GenMatrices.FaceMatrix;

    %% Initialization of the volume matrices for the trilinear components
    Matrices.Volume3L.M_NL_loc = GenMatrices.VolMatrix3L;

end
    
