%> @file   MatrixPreallocationPoroelasticity.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the matrices structure for the poroelasticity problem.
%>
%==========================================================================
%> @section classMatrixPreallocationPoroelasticity Class description
%==========================================================================
%> @brief           Preallocation of the matrices structure for poroelasticity 
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

function [Matrices] = MatrixPreallocationPoroelasticity(GenMatrices)

    %% Initialization of the volume matrices
    Matrices.Volume.V1_loc = GenMatrices.VolMatrix;
    Matrices.Volume.V2_loc = GenMatrices.VolMatrix;
    Matrices.Volume.V3_loc = GenMatrices.VolMatrix;
    Matrices.Volume.V4_loc = GenMatrices.VolMatrix;

    Matrices.Volume.M1_P_rho_loc  = GenMatrices.VolMatrix;
    Matrices.Volume.M1_P_rhof_loc = GenMatrices.VolMatrix;
    Matrices.Volume.M1_P_rhow_loc = GenMatrices.VolMatrix;

    Matrices.Volume.M1_P_eta_kper_loc  = GenMatrices.VolMatrix;
    Matrices.Volume.M1_P_rho2_zeta_loc = GenMatrices.VolMatrix;
    Matrices.Volume.M1_P_rho_zeta2_loc = GenMatrices.VolMatrix;

    Matrices.Volume.MPrjP_1_loc  = GenMatrices.VolMatrix;
  
    Matrices.Volume.B1_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B2_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B3_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B4_loc = GenMatrices.VolMatrix;

    Matrices.Volume.B1_beta_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B2_beta_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B3_beta_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B4_beta_loc = GenMatrices.VolMatrix;

    Matrices.Volume.B1_beta2_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B2_beta2_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B3_beta2_loc = GenMatrices.VolMatrix;
    Matrices.Volume.B4_beta2_loc = GenMatrices.VolMatrix;

    Matrices.Volume.P1_loc = GenMatrices.VolMatrix;
    Matrices.Volume.P2_loc = GenMatrices.VolMatrix;

    Matrices.Volume.P1_beta_loc = GenMatrices.VolMatrix;
    Matrices.Volume.P2_beta_loc = GenMatrices.VolMatrix;
    
    %% Initialization of the face matrices
    Matrices.Faces.IT1_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.IT2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.IT3_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.IT4_loc = GenMatrices.FaceMatrix;
        
    Matrices.Faces.S1_P_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S4_P_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.S1_B_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S2_B_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S3_B_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S4_B_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.S1_B_beta_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S2_B_beta_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S3_B_beta_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S4_B_beta_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.S1_B_beta2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S2_B_beta2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S3_B_beta2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.S4_B_beta2_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.BT1_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT3_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT4_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.BT1_beta_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT2_beta_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT3_beta_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT4_beta_loc = GenMatrices.FaceMatrix;

    Matrices.Faces.BT1_beta2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT2_beta2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT3_beta2_loc = GenMatrices.FaceMatrix;
    Matrices.Faces.BT4_beta2_loc = GenMatrices.FaceMatrix;

    %% Initialization of the absorbing boundary matrices
    Matrices.Faces.ABC_uu_1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_uu_2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_uu_3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_uu_4_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.ABC_uw_1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_uw_2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_uw_3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_uw_4_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.ABC_wu_1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_wu_2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_wu_3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_wu_4_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.ABC_ww_1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_ww_2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_ww_3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.ABC_ww_4_loc = GenMatrices.FaceMatrixBd;

    %% Initialization of the acoustic boundary matrices
    Matrices.Faces.D1_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.D2_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.D3_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.D4_loc = GenMatrices.FaceMatrixBd;

    %% Initialization of the elastic boundary matrices
    Matrices.Faces.S1_P_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S4_P_EP_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.S1_B_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S2_B_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S3_B_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S4_B_EP_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.S1_B_beta_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S2_B_beta_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S3_B_beta_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S4_B_beta_EP_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.S1_B_beta2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S2_B_beta2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S3_B_beta2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.S4_B_beta2_EP_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.BT1_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT3_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT4_EP_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.BT1_beta_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT2_beta_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT3_beta_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT4_beta_EP_loc = GenMatrices.FaceMatrixBd;

    Matrices.Faces.BT1_beta2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT2_beta2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT3_beta2_EP_loc = GenMatrices.FaceMatrixBd;
    Matrices.Faces.BT4_beta2_EP_loc = GenMatrices.FaceMatrixBd;

end
    