%> @file   VolumeMatricesAssemblyHeatQF.m
%> @author Mattia Corti
%> @date   5 February 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> heat equation (quadrature-free implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyHeatQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for heat equation (quadrature-free implementation).
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

function [Matrices] = VolumeMatricesAssemblyHeatQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        mu    = Data.mu{id}(0,0);
        sigma = Data.sigma{id}(0,0);
        
        % Mass matrices assembly
        Matrices.Mprj_loc(1:nbases^2,1) = Integral.phiphiC;
    
        Matrices.M_loc(1:nbases^2,1) = sigma*Integral.phiphiC;
        
        % Stiffness matrix assembly
        Matrices.A_loc(1:nbases^2,1) = mu*(Integral.gradxgradxC + Integral.gradygradyC);
 
end
                    