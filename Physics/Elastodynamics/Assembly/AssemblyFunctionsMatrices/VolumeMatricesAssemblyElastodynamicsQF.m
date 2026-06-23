%> @file   VolumeMatricesAssemblyElastodynamicsQF.m
%> @author Mattia Corti
%> @date   8 May 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> elastodynamics problem (quadrature-free implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyElastodynamicsQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for elastodynamics problem (quadrature-free implementation).
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

function [Matrices] = VolumeMatricesAssemblyElastodynamicsQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        rho_e = Data.rho_el{id}(0,0);
        zeta  = Data.zeta{id}(0,0);
        mu    = Data.mu_el{id}(0,0);
        lam   = Data.lam_el{id}(0,0);
            
        % Stiffness matrix assembly
        Matrices.V1_loc(1:nbases^2) = (lam+2*mu)*Integral.gradxgradxC + mu*Integral.gradygradyC;
        Matrices.V2_loc(1:nbases^2) = lam*Integral.gradygradxC + mu*Integral.gradxgradyC;
        Matrices.V3_loc(1:nbases^2) = lam*Integral.gradxgradyC + mu*Integral.gradygradxC;
        Matrices.V4_loc(1:nbases^2) = (lam+2*mu)*Integral.gradygradyC + mu*Integral.gradxgradxC;

        % Mass matrices assembly
        Matrices.M1_P_rho_loc(1:nbases^2) = rho_e*Integral.phiphiC;
        Matrices.MPrjP_1_loc(1:nbases^2)  = Integral.phiphiC;

        % Newmark matrices assembly
        Matrices.D1_loc(1:nbases^2)  = (2*rho_e*zeta)*Integral.phiphiC;
        Matrices.C1_loc(1:nbases^2)  = (rho_e*zeta^2)*Integral.phiphiC;
 
end
                    