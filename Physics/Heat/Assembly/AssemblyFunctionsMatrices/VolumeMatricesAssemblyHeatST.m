%> @file   VolumeMatricesAssemblyHeatST.m
%> @author Mattia Corti
%> @date   5 February 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> heat equation (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyHeatST description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for heat equation (subtriangulation implementation).
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

function [Matrices] = VolumeMatricesAssemblyHeatST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        mu    = Data.mu{id}(elem.xq,elem.yq);
        sigma = Data.sigma{id}(elem.xq,elem.yq);

        % Mass matrices assembly
        Matrices.Mprj_loc(1:nbases,1:nbases)  = (elem.dx .* elem.phiq)' * elem.phiq;
        
        Matrices.M_loc(1:nbases,1:nbases)  = (elem.dx .* sigma .* elem.phiq)' * elem.phiq;
        
        % Stiffness matrix assembly
        Matrices.A_loc(1:nbases,1:nbases)     = (elem.dx .* (mu .* elem.gradqx))' * elem.gradqx ...
                                              + (elem.dx .* (mu .* elem.gradqy))' * elem.gradqy;

end
                    
