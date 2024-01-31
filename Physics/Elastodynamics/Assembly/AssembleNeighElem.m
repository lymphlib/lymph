%> @file AssembleNeighElem.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Assembles the local matrices corresponding
%> to the neighbors of a given element.
%>
%==========================================================================
%> @section classAssembleNeighElem Class description
%==========================================================================
%> @brief            Assembles the local matrices for neighbors
%
%> @param M          Matrix to be assembled
%> @param row        Indeces of the current rows
%> @param M1         Local matrix of the neighboring element
%> @param nln        Number of basis functions
%> @param nE         Number of elements
%
%> @retval M         Matrix assembles
%>
%==========================================================================

function [M] = AssembleNeighElem(M, row, neight, M1, nln, nS, nE)

% % app = (neight > -1) .* (neight <= nE);
% % neight(app==0) = [];
% % j = ones(nln,1)*(neight-1)*nln + [1:nln]';
% % M1(:,:,app==0) = [];
% % for iedg = 1:size(M1,3)
% %     M(row,j(:,iedg)) = M(row,j(:,iedg)) + M1(:,:,iedg);
% % end
% app = (neight > -1) .* (neight >= nS + 1) .* (neight <= nE);
% neight(app==0) = [];
% % j = (ones(nln,1)*(neight-1)*nln + [1:nln]') - nS;
% j = (ones(nln,1)*(neight-nS-1)*nln + (1:nln)');
% M1(:,:,app==0) = [];
% for iedg = 1:size(M1,3)
%     M(row,j(:,iedg)) = M(row,j(:,iedg)) + M1(:,:,iedg);
% end

app = (neight > -1) .* (neight >= nS + 1) .* (neight <= nE);
neight(app==0) = [];
j = ones(nln,1)*(neight-1)*nln + [1:nln]';
M1(:,:,app==0) = [];
M1 = reshape(M1, size(M1,1), size(M1,2)*size(M1,3));
M(row,j) = M(row,j) + M1;
