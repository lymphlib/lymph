%> @file   VolumeForcingAssemblyHeatST.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   5 February 2026
%> @brief Assembly of the vectors associated with volume integrals for
%> heat equation (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyHeatST Class description
%==========================================================================
%> @brief            Assembly of the vectors associated with volume
%> integrals for heat equation (subtriangulation implementation).
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

function [Forcing] = VolumeForcingAssemblyHeatST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        f = Data.source_f{id}(elem.xq,elem.yq,AssembInfo.t);

        % Forcing term assembly
        Forcing.F_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * f;
        
end

