%> @file  AllocateMatrixElaLocalEl.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Allocation of local volume matrices 
%>
%==========================================================================
%> @section classAllocateMatrixElaLocalEl Class description
%==========================================================================
%> @brief Allocation of local volume matrices for the elastic domain 
%
%> @param nbases Number of basis function
%
%> @retval El    Struct with local matrices 
%>
%==========================================================================

function  [El] = AllocateMatrixElaLocalEl(nbases)

El.V1_loc       = zeros(nbases, nbases);
El.V2_loc       = zeros(nbases, nbases);
El.V3_loc       = zeros(nbases, nbases);
El.V4_loc       = zeros(nbases, nbases);
El.D1_loc       = zeros(nbases, nbases);
El.C1_loc       = zeros(nbases, nbases);
El.M1_P_rho_loc = zeros(nbases, nbases);
El.MPrjP_1_loc  = zeros(nbases, nbases);
