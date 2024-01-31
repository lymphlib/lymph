%> @file  LocalToGlobalMatrixElaEdgeNeigh.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Assembly of the edge matrices for the elastic domain
%(neighbouring elements)
%>
%==========================================================================
%> @section classLocalToGlobalMatrixElaEdgeNeigh Class description
%==========================================================================
%> @brief Assembly of the edge matrices for the elastic domain
%(neighbouring elements)
%
%> @param A         Struct with global matrices
%> @param Edge      Struct with local edge matrices
%> @param index     Vector with local to global mapping
%> @param neigh_ie  Vector with neigh element connectivity
%> @param nbases    number of local basis functions
%> @param nel_s     starting index
%> @param nel_e     ending index 
%
%> @retval A   Struct with global matrices
%>
%==========================================================================

function [A] = LocalToGlobalMatrixElaEdgeNeigh(A, Edge, index, neigh_ie, nbases, nel_s, nel_e)

% Assembling boundary DG matrices (Interior penalty and stabilization)
[A.S1_P]       = AssembleNeighElem(A.S1_P, index, neigh_ie, Edge.SN1_loc, nbases, nel_s, nel_e);
[A.S4_P]       = AssembleNeighElem(A.S4_P, index, neigh_ie, Edge.SN4_loc, nbases, nel_s, nel_e);

[A.IT1_P]      = AssembleNeighElem(A.IT1_P, index, neigh_ie, Edge.ITN1_loc, nbases, nel_s, nel_e);
[A.IT2_P]      = AssembleNeighElem(A.IT2_P, index, neigh_ie, Edge.ITN2_loc, nbases, nel_s, nel_e);
[A.IT3_P]      = AssembleNeighElem(A.IT3_P, index, neigh_ie, Edge.ITN3_loc, nbases, nel_s, nel_e);
[A.IT4_P]      = AssembleNeighElem(A.IT4_P, index, neigh_ie, Edge.ITN4_loc, nbases, nel_s, nel_e);



