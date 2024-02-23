%> @file  AllocateMatrixFluidGlobal.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Allocation of global matrices 
%>
%==========================================================================
%> @section classAllocateMatrixFluidGlobal Class description
%==========================================================================
%> @brief Allocation of global matrices for the fluid domain 
%>
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%>
%> @retval A   Struct with Global matrices
%>
%==========================================================================

function [A] = AllocateMatrixFluidGlobal(femregion, neighbor) 

% Mass Matrices
% \int_{\Omega_f}  mu^(-1) . (dev(sigma) . dev(tau) ) dx}
A.M_FF  = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_f*femregion.nbases^2);
A.MF_1 = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_f*femregion.nbases^2);
A.MF_2 = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_f*femregion.nbases^2);

% Projection Matrix
% \int_{\Omega_f} (sigma . tau ) dx
A.Mprj = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_f*femregion.nbases^2);

% Null Matrix
A.Z = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_f*femregion.nbases^2);

% Stiffness Matrices
% \int_{\Omega_f}   div(sigma) . div(tau) dx
A.B1  = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_p*femregion.nbases^2);
A.B2  = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_p*femregion.nbases^2);
A.B3  = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_p*femregion.nbases^2);
A.B4  = spalloc(femregion.ndof_f, femregion.ndof_f,femregion.nel_p*femregion.nbases^2);

% \int_{E_hp}   {div sigma} . [tau . n] ds
A.BT1  = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.BT2 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.BT3 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.BT4 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));

% \int_{E_hp} penalty_E  h_s^(-1) [sigma n].[tau n] ds  
A.S1_B = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.S2_B = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.S3_B = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.S4_B = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));

% \int_{Gamma_I}   (sigma n . n)(tau n . n) ds
A.C1 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C2 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C3 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C4 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C5 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C6 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C7 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C8 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
A.C9 = spalloc(femregion.ndof_f,femregion.ndof_f,femregion.nbases^2*sum(neighbor.nedges+1));
