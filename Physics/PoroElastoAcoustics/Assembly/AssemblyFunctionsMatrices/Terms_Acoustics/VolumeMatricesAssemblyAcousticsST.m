%> @file   VolumeMatricesAssemblyAcousticsST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   5 June 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> acoustics problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyAcousticsST description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for acoustics problem (subtriangulation implementation).
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

function [Matrices] = VolumeMatricesAssemblyAcousticsST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        rho_a   = Data.rho_a{id}(elem.xq,elem.yq);
        c       = Data.c{id}(elem.xq,elem.yq);
        
        % Matrices assembly
        Matrices.W_loc(1:nbases,1:nbases) = (elem.dx .* (rho_a .* elem.gradqx))' * elem.gradqx ...
                                          + (elem.dx .* (rho_a .* elem.gradqy))' * elem.gradqy;

        Matrices.M_P_A_loc(1:nbases,1:nbases) = (elem.dx .* (c.^(-2).* rho_a  .* elem.phiq))' * elem.phiq;
        Matrices.MPrjA_loc(1:nbases,1:nbases) = (elem.dx .* elem.phiq)' * elem.phiq;

end
                    
