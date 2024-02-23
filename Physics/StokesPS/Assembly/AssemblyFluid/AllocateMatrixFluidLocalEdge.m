%> @file  AllocateMatrixFluidLocalEdge.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Allocation of local edge matrices 
%>
%==========================================================================
%> @section classAllocateMatrixFluidLocalEdge Class description
%==========================================================================
%> @brief Allocation of local volume matrices for the fluid domain 
%>
%> @param nbases Number of basis function
%> @param nedges_ie_unq Number of edges of the polygonal element
%>
%> @retval Edge  Struct with local edge matrices 
%>
%==========================================================================

function [Edge] = AllocateMatrixFluidLocalEdge(nbases, nedges_ie_unq)

% Local matrix allocation for element itself
Edge.S1_B_loc          = zeros(nbases, nbases);
Edge.S2_B_loc          = zeros(nbases, nbases);
Edge.S3_B_loc          = zeros(nbases, nbases);
Edge.S4_B_loc          = zeros(nbases, nbases);
Edge.BT1_loc           = zeros(nbases, nbases);
Edge.BT2_loc           = zeros(nbases, nbases);
Edge.BT3_loc           = zeros(nbases, nbases);
Edge.BT4_loc           = zeros(nbases, nbases);
Edge.C1_loc            = zeros(nbases, nbases);
Edge.C2_loc            = zeros(nbases, nbases);
Edge.C3_loc            = zeros(nbases, nbases);
Edge.C4_loc            = zeros(nbases, nbases);
Edge.C5_loc            = zeros(nbases, nbases);
Edge.C6_loc            = zeros(nbases, nbases);
Edge.C7_loc            = zeros(nbases, nbases);
Edge.C8_loc            = zeros(nbases, nbases);
Edge.C9_loc            = zeros(nbases, nbases);

% Local matrix allocation for neighbors
% Edge.SN_B1_loc       = zeros(nbases, nbases, nedges_ie);
% Edge.SN_B2_loc       = zeros(nbases, nbases, nedges_ie);
% Edge.SN_B3_loc       = zeros(nbases, nbases, nedges_ie);
% Edge.SN_B4_loc       = zeros(nbases, nbases, nedges_ie);
% Edge.BTN1_loc        = zeros(nbases, nbases, nedges_ie);
% Edge.BTN2_loc        = zeros(nbases, nbases, nedges_ie);
% Edge.BTN3_loc        = zeros(nbases, nbases, nedges_ie);
% Edge.BTN4_loc        = zeros(nbases, nbases, nedges_ie);

Edge.SN_B1_loc       = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.SN_B2_loc       = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.SN_B3_loc       = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.SN_B4_loc       = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.BTN1_loc        = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.BTN2_loc        = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.BTN3_loc        = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.BTN4_loc        = zeros(nbases, nbases, length(nedges_ie_unq));
