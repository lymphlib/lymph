%> @file  AllocateMatrixElaLocalEdge.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Allocation of local edge matrices 
%>
%==========================================================================
%> @section classAllocateMatrixElaLocalEdge Class description
%==========================================================================
%> @brief Allocation of local volume matrices for the elastic domain 
%
%> @param nbases Number of basis function
%> @param nedges_ie_unq Number of edges of the polygonal element
%
%> @retval Edge  Struct with local edge matrices 
%>
%==========================================================================

function [Edge] = AllocateMatrixElaLocalEdge(nbases, nedges_ie_unq)

% Local matrix allocation for element itself
Edge.S1_P_loc  = zeros(nbases, nbases);
Edge.S4_P_loc  = zeros(nbases, nbases);

Edge.IT1_P_loc = zeros(nbases, nbases);
Edge.IT2_P_loc = zeros(nbases, nbases);
Edge.IT3_P_loc = zeros(nbases, nbases);
Edge.IT4_P_loc = zeros(nbases, nbases);

% Absorbing boundary
Edge.ABC_R1_loc = zeros(nbases, nbases);
Edge.ABC_R2_loc = zeros(nbases, nbases);
Edge.ABC_R3_loc = zeros(nbases, nbases);
Edge.ABC_R4_loc = zeros(nbases, nbases);

Edge.ABC_S1_loc = zeros(nbases, nbases);
Edge.ABC_S2_loc = zeros(nbases, nbases);
Edge.ABC_S3_loc = zeros(nbases, nbases);
Edge.ABC_S4_loc = zeros(nbases, nbases);

% % Local matrix allocation for neighbors
% Edge.ITN1_loc = zeros(nbases, nbases, nedges_ie);
% Edge.ITN2_loc = zeros(nbases, nbases, nedges_ie);
% Edge.ITN3_loc = zeros(nbases, nbases, nedges_ie);
% Edge.ITN4_loc = zeros(nbases, nbases, nedges_ie);
% 
% Edge.SN1_loc  = zeros(nbases, nbases, nedges_ie);
% Edge.SN4_loc  = zeros(nbases, nbases, nedges_ie);

Edge.ITN1_loc = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.ITN2_loc = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.ITN3_loc = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.ITN4_loc = zeros(nbases, nbases, length(nedges_ie_unq));

Edge.SN1_loc  = zeros(nbases, nbases, length(nedges_ie_unq));
Edge.SN4_loc  = zeros(nbases, nbases, length(nedges_ie_unq));
