%> @file  ReshapeVectorsMatrices.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 29 May 2026
%> @brief The function reshapes the vectors or the matrices accordingly to
%> the desired output form.
%>
%==========================================================================
%> @section classReshapeVectorsMatrices Class description
%==========================================================================
%> @brief            The function reshapes the vectors or the matrices accordingly to
%> the desired output form.
%>
%> @param M         Initial matrix or vector from the Assembly function
%> @param ii_index  Vector containing the first position in the final matrix
%> @param jj_index  Vector containing the second position in the final matrix
%> @param ii_index_neigh Vector containing the first position in the final matrix for faces matrices
%> @param jj_index_neigh Vector containing the second position in the final matrix for faces matrices
%> @param femregion Finite Element structure (see CreateDOF.m)
%> @param select    Contains the rows to maintain in the matrix assembly 
%> (typically the allocation is overestimated)
%> @param select_neigh Contains the rows to maintain in the matrix assembly
%> (typically the allocation is overestimated)
%> @param selectVec Contains the rows to maintain in the vector assembly
%> (typically the allocation is overestimated)
%> @param MatTag    Tag associated with the physics in the structure for multiphysics matrix assembly
%>
%> @retval M        Reshaped matrix or vector
%>
%==========================================================================

function [M] = ReshapeVectorsMatrices(M, ii_index, jj_index, ii_index_neigh, jj_index_neigh, femregion, select, select_neigh, selectVec, MatTag)

    M = vertcat(M{:});
    M = reshape(M, [1, size(M,1)*size(M,2)]);
    
    if length(M) == femregion.nel*max(femregion.nbases)
    
        % Vector assembly: delete the extra-spaces in matrices
        M = M(selectVec);
    
        % Vector assembly: vector transposition to column
        M = M';
    
    elseif length(M) == femregion.nel
    
        % Cell Vector assembly: vector transposition to column
        M = M';
    
    elseif length(M) == length(select)
    
        % Matrix assembly: delete the extra-spaces in matrices
        M = M(select);
    
        % Matrix assembly: sparse matrices construction
        M = sparse(ii_index, jj_index, M, femregion.ndof, femregion.ndof);
    
    else
    
        % Matrix assembly: delete the extra-spaces in matrices with neighbors elements
        M = M(select_neigh);
    
        % Matrix assembly: sparse matrices construction
        M = sparse(ii_index_neigh, jj_index_neigh, M, femregion.ndof, femregion.ndof);
    end

    if not(isempty(MatTag))
        idxmin_1 = sum(femregion.ndof_phys(1:find(femregion.phys == MatTag(1))-1))+1;
        idxmax_1 = sum(femregion.ndof_phys(1:find(femregion.phys == MatTag(1))));
        
        if size(M,2) == 1
            idxmin_2 = 1;
            idxmax_2 = 1;
        else
            idxmin_2 = sum(femregion.ndof_phys(1:find(femregion.phys == MatTag(2))-1))+1;
            idxmax_2 = sum(femregion.ndof_phys(1:find(femregion.phys == MatTag(2))));
        end

        M = M(idxmin_1:idxmax_1,idxmin_2:idxmax_2);
    end

end
