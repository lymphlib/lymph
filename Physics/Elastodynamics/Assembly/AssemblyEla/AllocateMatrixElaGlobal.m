%> @file  AllocateMatrixElaGlobal.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Allocation of global matrices 
%>
%==========================================================================
%> @section classAllocateMatrixElaGlobal Class description
%==========================================================================
%> @brief Allocation of global matrices for the elastic domain 
%
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%
%> @retval A   Struct with Global matrices
%>
%==========================================================================

function [A] = AllocateMatrixElaGlobal(femregion, neighbor) 

% Mass matrices
% \int_{\Omega_e}  rho_el (u . v ) dx}
A.M1_P_rho = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.M2_P_rho = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);

% Projection matrices
% \int_{\Omega_e}  (phi . psi ) dx} 
A.MPrjP_1 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.MPrjP_2 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);

% Damping matrices
% \int_{\Omega_e} (2 . rho . zeta  . u . v ) dx}
A.D1 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.D2 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
% \int_{\Omega_e} (rho . zeta^2 . u . v ) dx}
A.C1 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.C2 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);

% Absorbing matrices
A.ABC_S1 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.ABC_S2 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.ABC_S3 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.ABC_S4 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);

A.ABC_R1 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.ABC_R2 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.ABC_R3 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.ABC_R4 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);

% Stiffness matrices
% \int_{\Omega_e} (sigma(u) eps(v) dx}
A.V1 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.V2 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.V3 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.V4 = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
% \int_{E_he} penalty_E  h_s^(-1) [v].[u] ds
A.S1_P = spalloc(femregion.ndof_e,femregion.ndof_e,femregion.nbases^2*sum(neighbor.nedges+1));
A.S2_P = spalloc(femregion.ndof_e,femregion.ndof_e,femregion.nbases^2*sum(neighbor.nedges+1));
A.S3_P = spalloc(femregion.ndof_e,femregion.ndof_e,femregion.nbases^2*sum(neighbor.nedges+1));
A.S4_P = spalloc(femregion.ndof_e,femregion.ndof_e,femregion.nbases^2*sum(neighbor.nedges+1));
% \int_{E_he} {sigma(v)} . [u]ds
A.IT1_P = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.IT2_P = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.IT3_P = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);
A.IT4_P = spalloc(femregion.ndof_e, femregion.ndof_e,femregion.nel_e*femregion.nbases^2);

