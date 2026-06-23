%> @file   SIPMatricesPoroelasticity.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   12 May 2026
%> @brief  Global matrices structure for poroelasticity problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesPoroelasticity Class description
%==========================================================================
%> @brief           Global matrices structure for poroelasticity problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesPoroelasticity(Matrices_loc)

    %% Mass matrix elastic
    Matrices.M_P_rho  = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_rho_loc,Matrices_loc.Volume.M1_P_rho_loc);
    Matrices.M_P_rhof = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_rhof_loc,Matrices_loc.Volume.M1_P_rhof_loc);
    Matrices.M_P_rhow = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_rhow_loc,Matrices_loc.Volume.M1_P_rhow_loc);

    Matrices.M_P_eta_kper = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_eta_kper_loc,Matrices_loc.Volume.M1_P_eta_kper_loc);

    %% Projection matrix
    Matrices.MPrjP   = matlab.internal.math.blkdiag(Matrices_loc.Volume.MPrjP_1_loc,Matrices_loc.Volume.MPrjP_1_loc);
    Matrices.MPrjP_1 = Matrices_loc.Volume.MPrjP_1_loc;

    %% Damping matrix velocity
    Matrices.D_vel = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_rho2_zeta_loc,Matrices_loc.Volume.M1_P_rho2_zeta_loc);
    Matrices.D_dis = matlab.internal.math.blkdiag(Matrices_loc.Volume.M1_P_rho_zeta2_loc,Matrices_loc.Volume.M1_P_rho_zeta2_loc);

    %% Absorbing matrices
    Matrices.ABC_UU = [Matrices_loc.Faces.ABC_uu_1_loc, Matrices_loc.Faces.ABC_uu_2_loc; 
                       Matrices_loc.Faces.ABC_uu_3_loc, Matrices_loc.Faces.ABC_uu_4_loc];

    Matrices.ABC_UW = [Matrices_loc.Faces.ABC_uw_1_loc, Matrices_loc.Faces.ABC_uw_2_loc;
                       Matrices_loc.Faces.ABC_uw_3_loc, Matrices_loc.Faces.ABC_uw_4_loc];

    Matrices.ABC_WU = [Matrices_loc.Faces.ABC_wu_1_loc, Matrices_loc.Faces.ABC_wu_2_loc; 
                       Matrices_loc.Faces.ABC_wu_3_loc, Matrices_loc.Faces.ABC_wu_4_loc];

    Matrices.ABC_WW = [Matrices_loc.Faces.ABC_ww_1_loc, Matrices_loc.Faces.ABC_ww_2_loc; 
                       Matrices_loc.Faces.ABC_ww_3_loc, Matrices_loc.Faces.ABC_ww_4_loc];

    %% Imperfect pore matrix
    Matrices.D_P  = [Matrices_loc.Faces.D1_loc,    Matrices_loc.Faces.D2_loc;    
                     Matrices_loc.Faces.D3_loc,    Matrices_loc.Faces.D4_loc];
     
    %% Elastic SIP-DG matrix
    V    = [Matrices_loc.Volume.V1_loc,    Matrices_loc.Volume.V2_loc;    
            Matrices_loc.Volume.V3_loc,    Matrices_loc.Volume.V4_loc];
     
    IT_P = [Matrices_loc.Faces.IT1_loc, Matrices_loc.Faces.IT2_loc; 
            Matrices_loc.Faces.IT3_loc, Matrices_loc.Faces.IT4_loc];
    
    S_P  = matlab.internal.math.blkdiag(Matrices_loc.Faces.S1_P_loc,Matrices_loc.Faces.S4_P_loc);
    
    S_P_EP  = matlab.internal.math.blkdiag(Matrices_loc.Faces.S1_P_EP_loc,Matrices_loc.Faces.S4_P_EP_loc);
    
    Matrices.A_E = V + S_P + S_P_EP - IT_P - IT_P';
    Matrices.DGe = V + S_P;


    B    = [Matrices_loc.Volume.B1_loc,    Matrices_loc.Volume.B2_loc;    
            Matrices_loc.Volume.B3_loc,    Matrices_loc.Volume.B4_loc];
    
    BT = [Matrices_loc.Faces.BT1_loc, Matrices_loc.Faces.BT2_loc; 
          Matrices_loc.Faces.BT3_loc, Matrices_loc.Faces.BT4_loc];
    
    S_B  = [Matrices_loc.Faces.S1_B_loc, Matrices_loc.Faces.S2_B_loc; 
            Matrices_loc.Faces.S3_B_loc, Matrices_loc.Faces.S4_B_loc];

    BT_EP = [Matrices_loc.Faces.BT1_EP_loc, Matrices_loc.Faces.BT2_EP_loc;
             Matrices_loc.Faces.BT3_EP_loc, Matrices_loc.Faces.BT4_EP_loc];
    
    S_B_EP = [Matrices_loc.Faces.S1_B_EP_loc, Matrices_loc.Faces.S2_B_EP_loc; 
              Matrices_loc.Faces.S3_B_EP_loc, Matrices_loc.Faces.S4_B_EP_loc];

    Matrices.A_P         = B - BT - BT' + S_B - BT_EP - BT_EP' + S_B_EP;
    Matrices.DGp = B + S_B;


    B_beta = [Matrices_loc.Volume.B1_beta_loc,    Matrices_loc.Volume.B2_beta_loc;    
              Matrices_loc.Volume.B3_beta_loc,    Matrices_loc.Volume.B4_beta_loc];
    
    BT_beta = [Matrices_loc.Faces.BT1_beta_loc, Matrices_loc.Faces.BT2_beta_loc; 
               Matrices_loc.Faces.BT3_beta_loc, Matrices_loc.Faces.BT4_beta_loc];
    
    S_B_beta  = [Matrices_loc.Faces.S1_B_beta_loc, Matrices_loc.Faces.S2_B_beta_loc; 
                 Matrices_loc.Faces.S3_B_beta_loc, Matrices_loc.Faces.S4_B_beta_loc];

    BT_beta_EP = [Matrices_loc.Faces.BT1_beta_EP_loc, Matrices_loc.Faces.BT2_beta_EP_loc;
                  Matrices_loc.Faces.BT3_beta_EP_loc, Matrices_loc.Faces.BT4_beta_EP_loc];
    
    S_B_beta_EP = [Matrices_loc.Faces.S1_B_beta_EP_loc, Matrices_loc.Faces.S2_B_beta_EP_loc; 
                   Matrices_loc.Faces.S3_B_beta_EP_loc, Matrices_loc.Faces.S4_B_beta_EP_loc];

    Matrices.A_P_beta_pf = B_beta - BT_beta - BT_beta' + S_B_beta - BT_beta_EP - BT_beta_EP' + S_B_beta_EP;
    Matrices.A_P_beta_fp = B_beta - BT_beta - BT_beta' + S_B_beta - BT_beta_EP - BT_beta_EP' + S_B_beta_EP;
    Matrices.DGp_beta = B_beta + S_B_beta;
   
    
    B_beta2 = [Matrices_loc.Volume.B1_beta2_loc,    Matrices_loc.Volume.B2_beta2_loc;    
               Matrices_loc.Volume.B3_beta2_loc,    Matrices_loc.Volume.B4_beta2_loc];
    
    BT_beta2 = [Matrices_loc.Faces.BT1_beta2_loc, Matrices_loc.Faces.BT2_beta2_loc; 
                Matrices_loc.Faces.BT3_beta2_loc, Matrices_loc.Faces.BT4_beta2_loc];
    
    S_B_beta2  = [Matrices_loc.Faces.S1_B_beta2_loc, Matrices_loc.Faces.S2_B_beta2_loc; 
                  Matrices_loc.Faces.S3_B_beta2_loc, Matrices_loc.Faces.S4_B_beta2_loc];

    BT_beta2_EP = [Matrices_loc.Faces.BT1_beta2_EP_loc, Matrices_loc.Faces.BT2_beta2_EP_loc;
                   Matrices_loc.Faces.BT3_beta2_EP_loc, Matrices_loc.Faces.BT4_beta2_EP_loc];
    
    S_B_beta2_EP = [Matrices_loc.Faces.S1_B_beta2_EP_loc, Matrices_loc.Faces.S2_B_beta2_EP_loc; 
                    Matrices_loc.Faces.S3_B_beta2_EP_loc, Matrices_loc.Faces.S4_B_beta2_EP_loc];
 
    Matrices.A_P_beta2 = B_beta2 - BT_beta2 - BT_beta2' + S_B_beta2 - BT_beta2_EP - BT_beta2_EP' + S_B_beta2_EP;
   
    %% Coupling matrices
    Matrices.P1 = Matrices_loc.Volume.P1_loc;
    Matrices.P2 = Matrices_loc.Volume.P2_loc;

    Matrices.P1_beta = Matrices_loc.Volume.P1_beta_loc;
    Matrices.P2_beta = Matrices_loc.Volume.P2_beta_loc;

end