%> @file   VolumeMatricesAssemblyElastodynamicsST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   8 May 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> elastodynamics problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyElastodynamicsST description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for elastodynamics problem (subtriangulation implementation).
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

function [Matrices] = VolumeMatricesAssemblyElastodynamicsST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        rho_e = Data.rho_el{id}(elem.xq,elem.yq);
        zeta  = Data.zeta{id}(elem.xq,elem.yq);
        mu    = Data.mu_el{id}(elem.xq,elem.yq);
        lam   = Data.lam_el{id}(elem.xq,elem.yq);
        
        % Mass matrices assembly
        Matrices.M1_P_rho_loc(1:nbases,1:nbases) =  (elem.dx .* (rho_e .* elem.phiq))' * elem.phiq;
        Matrices.MPrjP_1_loc(1:nbases,1:nbases)  = (elem.dx .* elem.phiq)' * elem.phiq;

        % Stiffness matrix assembly
        Matrices.V1_loc(1:nbases,1:nbases) = (elem.dx .* ((lam+2*mu) .* elem.gradqx))' * elem.gradqx ...
                                               + (elem.dx .* (mu  .* elem.gradqy))' * elem.gradqy;

        Matrices.V2_loc(1:nbases,1:nbases) = (elem.dx .* (lam .* elem.gradqx))' * elem.gradqy ...
                                               + (elem.dx .* (mu  .* elem.gradqy))' * elem.gradqx;

        Matrices.V3_loc(1:nbases,1:nbases) = (elem.dx .* (lam .* elem.gradqy))' * elem.gradqx ...
                                               + (elem.dx .* (mu  .* elem.gradqx))' * elem.gradqy;

        Matrices.V4_loc(1:nbases,1:nbases) = (elem.dx .* ((lam+2*mu) .* elem.gradqy))' * elem.gradqy ...
                                               + (elem.dx .* (mu  .* elem.gradqx))' * elem.gradqx;

        % Newmark matrices assembly
        Matrices.D1_loc(1:nbases,1:nbases) =  (elem.dx .* (2 * rho_e .* zeta .* elem.phiq))' * elem.phiq;
        Matrices.C1_loc(1:nbases,1:nbases) =  (elem.dx .* (2 * rho_e .* zeta.^2 .* elem.phiq))' * elem.phiq;
       
end
                    
