%> @file  AssembleODEMatricesPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 26 June 2024
%> @brief  Assembly of the matrices A,B and C for the problem Ax'' + Bx' + Cx = F
%>
%==========================================================================
%> @section classAssembleODEMatricesPoroAcuEla Class description
%> @brief  Assembly of the matrices A,B and C for the problem Ax'' + Bx' + Cx = F
%
%> @param Data       Struct with problem's data
%> @param Matrices   Matrices for poro-acoustic-elastic domain
%
%> @retval A,B,C     Matrices for the problem Ax'' + Bx' + Cx = F
%>
%==========================================================================

function [A, B, C] = AssembleODEMatricesPoroAcuEla(Data,Matrices)  
                            

[~,np1] = size(Matrices.Poro.M_P_rho);
[~,np2] = size(Matrices.Poro.M_P_rhow);
[~,nac] = size(Matrices.Acu.M_A);
[~,nel] = size(Matrices.Ela.M_P_rho);
 
[r,s] = size(Matrices.PoroAcu.C1_P);
[~,m] = size(Matrices.ElaAcu.C1_A);



A = [Matrices.Poro.M_P_rho,   Matrices.Poro.M_P_rhof,  sparse(np1,nac),   sparse(np1,nel);
     Matrices.Poro.M_P_rhof,  Matrices.Poro.M_P_rhow,  sparse(np2,nac),   sparse(np2,nel);
     sparse(nac,np1),         sparse(nac,np2),         Matrices.Acu.M_A,  sparse(nac,nel);
     sparse(nel,np1),         sparse(nel,np2),         sparse(nel,nac),   Matrices.Ela.M_P_rho;];
 

C = [Matrices.Poro.A_E + Matrices.Poro.A_P_beta2 + Matrices.Poro.D_dis, Matrices.Poro.A_P_beta_pf,    sparse(np1,nac),    Matrices.PoroEla.C_P;    
     Matrices.Poro.A_P_beta_fp,                                         Matrices.Poro.A_P,            sparse(np2,nac),    sparse(np2,nel); 
     sparse(nac,np1),                                                   sparse(nac,np2),              Matrices.Acu.A_A,   sparse(nac,nel);
     Matrices.PoroEla.C_E,                                              sparse(nel,np2),              sparse(nel,nac),    Matrices.Ela.A_E + Matrices.Ela.Ddis+Matrices.Ela.Rdis];

if (Data.tau ~=0)
    
  D = Matrices.Poro.M_P_eta_kper + (1-Data.tau)/Data.tau*Matrices.Poro.D_P;
  
  B = [Matrices.Poro.ABC_UU + Matrices.Poro.D_vel,  Matrices.Poro.ABC_UW,      Matrices.PoroAcu.C1_P,    sparse(r,m);
       Matrices.Poro.ABC_WU,                        Matrices.Poro.ABC_WW + D,  Matrices.PoroAcu.C2_P;    sparse(r,m);
       Matrices.PoroAcu.C1_A,                       Matrices.PoroAcu.C2_A,     sparse(s,s),              Matrices.ElaAcu.C1_A;
       sparse(m,r),                                 sparse(m,r)                Matrices.ElaAcu.C1_E      Matrices.Ela.Dvel + Matrices.Ela.Svel];
else
    
  D = Matrices.Poro.M_P_eta_kper;
    
  B = [Matrices.Poro.ABC_UU + Matrices.Poro.D_vel,  Matrices.Poro.ABC_UW,       Matrices.PoroAcu.C1_P,  sparse(r,m);
       Matrices.Poro.ABC_WU,                        Matrices.Poro.ABC_WW + D,   sparse(r,s),            sparse(r,m);
       Matrices.PoroAcu.C1_A,                       sparse(s,r),                sparse(s,s),            Matrices.ElaAcu.C1_A;
       sparse(m,r),                                 sparse(m,r)                 Matrices.ElaAcu.C1_E    Matrices.Ela.Dvel + Matrices.Ela.Svel];
        
 end


