%> @file  AllocateMatrixFluidLocalEl.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Allocation of local volume matrices 
%>
%==========================================================================
%> @section classAllocateMatrixFluidLocalEl Class description
%==========================================================================
%> @brief Allocation of local volume matrices for the poroelastic domain 
%>
%> @param nbases Number of basis function
%>
%> @retval El    Struct with local matrices 
%>
%==========================================================================

function [El] = AllocateMatrixFluidLocalEl(nbases)

El.MFF_loc      = zeros(nbases, nbases);
El.MF_1_loc     = zeros(nbases, nbases);
El.MF_2_loc     = zeros(nbases, nbases);
El.Mprj_loc     = zeros(nbases, nbases);
El.B1_loc       = zeros(nbases, nbases);
El.B2_loc       = zeros(nbases, nbases);
El.B3_loc       = zeros(nbases, nbases);
El.B4_loc       = zeros(nbases, nbases);
