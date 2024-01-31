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

function [M] = AssembleNeighEl(M, row, neigh_ie, MN, nbases, nStart, nEnd)

app = (neigh_ie > -1) .* (neigh_ie >= nStart) .* (neigh_ie <= nEnd);
neigh_ie(app==0) = [];
j = (ones(nbases,1)*(neigh_ie-nStart)*nbases + (1:nbases)');
MN(:,:,app==0) = [];
for iedg = 1:size(MN,3)
    M(row,j(:,iedg)) = M(row,j(:,iedg)) + MN(:,:,iedg);
end
