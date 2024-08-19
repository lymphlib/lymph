%> @file AssembleNeighEl.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Given a local matrix defined on an element, it adds the
%> corresponding local matrices defined on the element's neighbors to it.
%>
%==========================================================================
%> @section classAssembleNeighEl Class description
%==========================================================================
%> @brief            Assembles the local matrices for neighbors
%
%> @param M          Matrix to be assembled, starting from input [in, out]
%> @param row        Indices of the rows associated to current element
%> @param neigh_ie   Indices of the rows associated to neighboring elements
%> @param MN         Local matrices of the neighboring elements: third index
%>                   spans the neighboring elements by interface edge index
%> @param nbases     Number of basis functions
%> @param nStart     Starting index (typically = 1 in most single-physics systems)
%> @param nEnd       Ending index (typically = nStart + number of elements)
%>
%> @retval M         Matrix assembled
%>
%==========================================================================
% function [M] = AssembleNeighEl(M,row,neigh_ie,MN,nbases,nStart,nEnd)

function [M] = AssembleNeighEl(M, row, neight, MN, nbases, nStart, nEnd)

app = (neight > -1) .* (neight >= nStart) .* (neight <= nEnd);
neight(app==0) = [];
j = ones(nbases,1)*(neight-nStart)*nbases + (1:nbases)';
MN(:,:,app==0) = [];
MN = reshape(MN, size(MN,1), size(MN,2)*size(MN,3));
M(row,j) = M(row,j) + MN;
