%> @file   IPVolume3LMatricesAssemblyFKST.m
%> @author Mattia Corti
%> @date   15 September 2025
%> @brief Assembly of the matrices associated with volume integrals of trilinear forms for FK
%> problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classIPVolume3LMatricesAssemblyFKST Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals of trilinear forms for FK problem (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Matrices   Matrices struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = IPVolume3LMatricesAssemblyFKST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        alpha  = Data.alpha{id}(elem.xq,elem.yq,0);
                    
        % 3L mass matrix assembly
        M_NL_loc = zeros(nbases, nbases, nbases);

        for j = 1:nbases
            M_NL_loc(1:nbases,1:nbases,j) = ((elem.dx.* alpha .* elem.phiq(:,j)) .* elem.phiq )'*elem.phiq;
        end

        Matrices.M_NL_loc(1:nbases^2,1:nbases) = reshape(M_NL_loc,[nbases^2, nbases]);
end
                    
