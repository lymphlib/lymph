%> @file   VolumeMatricesAssemblyPoroelasticityST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   05 June 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> poroelasticity problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyPoroelasticityST description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for poroelasticity problem (subtriangulation implementation).
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

function [Matrices] = VolumeMatricesAssemblyPoroelasticityST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        mu        = Data.mu{id}(elem.xq,elem.yq);
        lam       = Data.lam{id}(elem.xq,elem.yq);
        beta      = Data.beta{id}(elem.xq,elem.yq);
        m         = Data.m{id}(elem.xq,elem.yq);
        rho_f     = Data.rho_f{id}(elem.xq,elem.yq);
        rho_w     = Data.rho_w{id}(elem.xq,elem.yq);
        rho       = Data.rho{id}(elem.xq,elem.yq);
        eta_kper  = Data.eta{id}(elem.xq,elem.yq)./Data.k_per{id}(elem.xq,elem.yq);
        rho2_zeta = 2*Data.rho{id}(elem.xq,elem.yq).*Data.zetap{id}(elem.xq,elem.yq);
        rho_zeta2 = Data.rho{id}(elem.xq,elem.yq).*Data.zetap{id}(elem.xq,elem.yq).^2;

        % Mass matrices assembly
        Matrices.MPrjP_1_loc(1:nbases,1:nbases)  = (elem.dx .* elem.phiq)' * elem.phiq;

        Matrices.M1_P_rhof_loc(1:nbases,1:nbases) =  (elem.dx .* (rho_f .* elem.phiq))' * elem.phiq;
        Matrices.M1_P_rhow_loc(1:nbases,1:nbases) =  (elem.dx .* (rho_w .* elem.phiq))' * elem.phiq;
        Matrices.M1_P_rho_loc(1:nbases,1:nbases)  =  (elem.dx .* (rho .* elem.phiq))' * elem.phiq;

        Matrices.M1_P_rho2_zeta_loc(1:nbases,1:nbases) =  (elem.dx .* (rho2_zeta .* elem.phiq))' * elem.phiq;
        Matrices.M1_P_rho_zeta2_loc(1:nbases,1:nbases) =  (elem.dx .* (rho_zeta2 .* elem.phiq))' * elem.phiq;
        Matrices.M1_P_eta_kper_loc(1:nbases,1:nbases)  =  (elem.dx .* (eta_kper .* elem.phiq))' * elem.phiq;

        % Stiffness matrix assembly
        Matrices.V1_loc(1:nbases,1:nbases) = (elem.dx .* ((lam+2*mu) .* elem.gradqx))' * elem.gradqx ...
                                           + (elem.dx .* (mu  .* elem.gradqy))' * elem.gradqy;

        Matrices.V2_loc(1:nbases,1:nbases) = (elem.dx .* (lam .* elem.gradqx))' * elem.gradqy ...
                                           + (elem.dx .* (mu  .* elem.gradqy))' * elem.gradqx;

        Matrices.V3_loc(1:nbases,1:nbases) = (elem.dx .* (lam .* elem.gradqy))' * elem.gradqx ...
                                           + (elem.dx .* (mu  .* elem.gradqx))' * elem.gradqy;

        Matrices.V4_loc(1:nbases,1:nbases) = (elem.dx .* ((lam+2*mu) .* elem.gradqy))' * elem.gradqy ...
                                           + (elem.dx .* (mu  .* elem.gradqx))' * elem.gradqx;

        Matrices.B1_loc(1:nbases,1:nbases) = (elem.dx .* (m .* elem.gradqx))' * elem.gradqx;
        Matrices.B2_loc(1:nbases,1:nbases) = (elem.dx .* (m .* elem.gradqx))' * elem.gradqy;
        Matrices.B3_loc(1:nbases,1:nbases) = (elem.dx .* (m .* elem.gradqy))' * elem.gradqx;
        Matrices.B4_loc(1:nbases,1:nbases) = (elem.dx .* (m .* elem.gradqy))' * elem.gradqy;

        Matrices.B1_beta_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta .* elem.gradqx))' * elem.gradqx;
        Matrices.B2_beta_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta .* elem.gradqx))' * elem.gradqy;
        Matrices.B3_beta_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta .* elem.gradqy))' * elem.gradqx;
        Matrices.B4_beta_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta .* elem.gradqy))' * elem.gradqy;

        Matrices.B1_beta2_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta.^2 .* elem.gradqx))' * elem.gradqx;
        Matrices.B2_beta2_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta.^2 .* elem.gradqx))' * elem.gradqy;
        Matrices.B3_beta2_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta.^2 .* elem.gradqy))' * elem.gradqx;
        Matrices.B4_beta2_loc(1:nbases,1:nbases) = (elem.dx .* (m .* beta.^2 .* elem.gradqy))' * elem.gradqy;

        Matrices.P1_loc(1:nbases,1:nbases) = - (elem.dx .* (m .* elem.phiq))' * elem.gradqx;
        Matrices.P2_loc(1:nbases,1:nbases) = - (elem.dx .* (m .* elem.phiq))' * elem.gradqy;

        Matrices.P1_beta_loc(1:nbases,1:nbases) = - (elem.dx .* (m .* beta .* elem.phiq))' * elem.gradqx;
        Matrices.P2_beta_loc(1:nbases,1:nbases) = - (elem.dx .* (m .* beta .* elem.phiq))' * elem.gradqy;
            
end
                    
