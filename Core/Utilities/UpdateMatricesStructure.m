%> @file  UpdateMatricesStructure.m
%> @author Mattia Corti
%> @date   29 May 2026
%> @brief  Conversion between element-wise and global matrices structures
%>
%==========================================================================
%> @section classUpdateMatricesStructure Class description
%==========================================================================
%> @brief Convert a cell array of element-wise matrices into the global
%>        Matrices struct format used in the assembly routines.
%>
%> @param Matrices_cell Cell array of length femregion.nel. Each entry
%>        Matrices_cell{ie} is a struct collecting all local matrices
%>        associated with element ie. 
%>
%> @retval Matrices  Struct in the global format. For each top-level field
%>        F, Matrices.F is a struct whose subfields are cell arrays indexed
%>        by the element.
%>
%==========================================================================

function Matrices = UpdateMatricesStructure(Matrices_cell)

    % Extraction of field names and preallocation
    matFields = fieldnames(Matrices_cell{1});
    Matrices = struct();

    %% Cycle over the struct fields
    for m = 1:numel(matFields)
        
        % Extract the subfield names and preallocation
        subFields = fieldnames(Matrices_cell{1}.(matFields{m}));
        SubStruct = struct();

        for k = 1:numel(subFields)
              
            if isstruct(Matrices_cell{1}.(matFields{m}).(subFields{k}))

                % Extract the subfield names and preallocation
                subsubFields = fieldnames(Matrices_cell{1}.(matFields{m}).(subFields{k}));

                for j = 1:numel(subsubFields)
                    % Preallocation of the k-th subfield
                    SubStruct.(subFields{k}).(subsubFields{j}) = cell(numel(Matrices_cell),1);

                    % Filling the k-th subfield
                    for ie = 1:numel(Matrices_cell)
                        SubStruct.(subFields{k}).(subsubFields{j}){ie} = Matrices_cell{ie}.(matFields{m}).(subFields{k}).(subsubFields{j});
                    end
                end
            else
                % Preallocation of the k-th subfield
                SubStruct.(subFields{k}) = cell(numel(Matrices_cell),1);
    
                % Filling the k-th subfield
                for ie = 1:numel(Matrices_cell)
                    SubStruct.(subFields{k}){ie} = Matrices_cell{ie}.(matFields{m}).(subFields{k});
                end
            end
        end

        % Final structure construction
        Matrices.(matFields{m}) = SubStruct;

    end
end
