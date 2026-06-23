%> @file   VolumeForcingAssemblyLaplacianST.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   12 May 2026
%> @brief Assembly of the vectors associated with volume integrals for
%> Laplacian problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyLaplacianST Class description
%==========================================================================
%> @brief            Assembly of the vectors associated with volume
%> integrals for Laplacian problem (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Forcing    Forcing struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = VolumeForcingAssemblyLaplacianST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        fSource = Data.source{id}(elem.xq,elem.yq);

        % Forcing term assembly
        Forcing.F_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * fSource;
        
end

