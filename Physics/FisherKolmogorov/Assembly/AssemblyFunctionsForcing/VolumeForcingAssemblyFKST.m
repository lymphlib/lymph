%> @file   VolumeForcingAssemblyFKST.m
%> @author Mattia Corti
%> @date   14 April 2026
%> @brief Assembly of the forcing terms with volume integrals for FK
%> problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyFKST Class description
%==========================================================================
%> @brief            Assembly of the forcing terms associated with volume 
%> integrals for FK problem (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Forcing    Matrices struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = VolumeForcingAssemblyFKST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        D_ext  = Data.D_ext{id}(elem.xq,elem.yq,AssembInfo.t);
        alpha  = Data.alpha{id}(elem.xq,elem.yq,AssembInfo.t);

        f      = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t,D_ext,alpha);

        % Forcing term assembly
        Forcing.F_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * f;
        
end
                    
