%> @file   VolumeMatricesAssemblyPoroelasticityQF.m
%> @author Mattia Corti
%> @date   29 May 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> poroelasticity problem (quadrature-free implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyPoroelasticityQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for poroelasticity problem (quadrature-free implementation).
%>
%> @param Data       Struct with problem's data
%> @param Integral   Integral struct containing values of bases integrals
%> @param Matrices   Matrices struct (to be modified)
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = VolumeMatricesAssemblyPoroelasticityQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        mu        = Data.mu{id}(0,0);
        lam       = Data.lam{id}(0,0);
        beta      = Data.beta{id}(0,0);
        m         = Data.m{id}(0,0);
        rho_f     = Data.rho_f{id}(0,0);
        rho_w     = Data.rho_w{id}(0,0);
        rho       = Data.rho{id}(0,0);
        eta_kper  = Data.eta{id}(0,0)./Data.k_per{id}(0,0);
        rho2_zeta = 2 * Data.rho{id}(0,0) .* Data.zetap{id}(0,0);
        rho_zeta2 = Data.rho{id}(0,0) .* Data.zetap{id}(0,0).^2;
            
        
        % Mass matrices assembly
        Matrices.M1_P_rhof_loc(1:nbases^2) = rho_f*Integral.phiphiC;
        Matrices.M1_P_rhow_loc(1:nbases^2) = rho_w*Integral.phiphiC;
        Matrices.M1_P_rho_loc(1:nbases^2)  = rho*Integral.phiphiC;
        
        Matrices.MPrjP_1_loc(1:nbases^2)   = Integral.phiphiC;

        Matrices.M1_P_rho2_zeta_loc(1:nbases^2) = rho2_zeta*Integral.phiphiC;
        Matrices.M1_P_rho_zeta2_loc(1:nbases^2) = rho_zeta2*Integral.phiphiC;
        Matrices.M1_P_eta_kper_loc(1:nbases^2)  = eta_kper*Integral.phiphiC;

        % Stiffness matrix assembly
        Matrices.V1_loc(1:nbases^2) = (lam+2*mu)*Integral.gradxgradxC + mu*Integral.gradygradyC;
        Matrices.V2_loc(1:nbases^2) = lam*Integral.gradygradxC + mu*Integral.gradxgradyC;
        Matrices.V3_loc(1:nbases^2) = lam*Integral.gradxgradyC + mu*Integral.gradygradxC;
        Matrices.V4_loc(1:nbases^2) = (lam+2*mu)*Integral.gradygradyC + mu*Integral.gradxgradxC;

        Matrices.B1_loc(1:nbases^2) = m*Integral.gradxgradxC;
        Matrices.B2_loc(1:nbases^2) = m*Integral.gradygradxC;
        Matrices.B3_loc(1:nbases^2) = m*Integral.gradxgradyC;
        Matrices.B4_loc(1:nbases^2) = m*Integral.gradygradyC;

        Matrices.B1_beta_loc(1:nbases^2) = (m .* beta)*Integral.gradxgradxC;
        Matrices.B2_beta_loc(1:nbases^2) = (m .* beta)*Integral.gradygradxC;
        Matrices.B3_beta_loc(1:nbases^2) = (m .* beta)*Integral.gradxgradyC;
        Matrices.B4_beta_loc(1:nbases^2) = (m .* beta)*Integral.gradygradyC;

        Matrices.B1_beta2_loc(1:nbases^2) = (m .* beta.^2)*Integral.gradxgradxC;
        Matrices.B2_beta2_loc(1:nbases^2) = (m .* beta.^2)*Integral.gradygradxC;
        Matrices.B3_beta2_loc(1:nbases^2) = (m .* beta.^2)*Integral.gradxgradyC;
        Matrices.B4_beta2_loc(1:nbases^2) = (m .* beta.^2)*Integral.gradygradyC;

        Matrices.P1_loc(1:nbases^2) = - m*Integral.gradxphiC;
        Matrices.P2_loc(1:nbases^2) = - m*Integral.gradyphiC;

        Matrices.P1_beta_loc(1:nbases^2) = - (m .* beta)*Integral.gradxphiC;
        Matrices.P2_beta_loc(1:nbases^2) = - (m .* beta)*Integral.gradyphiC;
 
end
                    