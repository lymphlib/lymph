%> @file  LocalToGlobalMatrixFluidEdgeNeigh.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Assembly of the edge matrices for the fluid domain
%(neighbouring elements)
%>
%==========================================================================
%> @section classLocalToGlobalMatrixFluidEdgeNeigh Class description
%==========================================================================
%> @brief Assembly of the edge matrices for the fluid domain
%(neighbouring elements)
%>
%> @param A         Struct with global matrices
%> @param Edge      Struct with local edge matrices
%> @param index     Vector with local to global mapping
%> @param neigh_ie  Vector with neigh element connectivity
%> @param nbases    number of local basis functions
%> @param nel_s     starting index
%> @param nel_e     ending index 
%> 
%> @retval A   Struct with global matrices
%>
%==========================================================================

function [A] = LocalToGlobalMatrixFluidEdgeNeigh(A, Edge, index, neigh_ie, nbases, nel_s, nel_e)

% Assembling boundary DG matrices (Interior penalty and stabilization)
[A.S1_B]       = AssembleNeighElem(A.S1_B,       index, neigh_ie, Edge.SN_B1_loc,       nbases, nel_s, nel_e);
[A.S2_B]       = AssembleNeighElem(A.S2_B,       index, neigh_ie, Edge.SN_B2_loc,       nbases, nel_s, nel_e);
[A.S3_B]       = AssembleNeighElem(A.S3_B,       index, neigh_ie, Edge.SN_B3_loc,       nbases, nel_s, nel_e);
[A.S4_B]       = AssembleNeighElem(A.S4_B,       index, neigh_ie, Edge.SN_B4_loc,       nbases, nel_s, nel_e);
[A.BT1]        = AssembleNeighElem(A.BT1,        index, neigh_ie, Edge.BTN1_loc,        nbases, nel_s, nel_e);
[A.BT2]        = AssembleNeighElem(A.BT2,        index, neigh_ie, Edge.BTN2_loc,        nbases, nel_s, nel_e);
[A.BT3]        = AssembleNeighElem(A.BT3,        index, neigh_ie, Edge.BTN3_loc,        nbases, nel_s, nel_e);
[A.BT4]        = AssembleNeighElem(A.BT4,        index, neigh_ie, Edge.BTN4_loc,        nbases, nel_s, nel_e);
