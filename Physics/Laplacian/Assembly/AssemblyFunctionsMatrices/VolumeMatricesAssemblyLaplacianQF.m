%> @file   VolumeMatricesAssemblyLaplacianQF.m
%> @author Mattia Corti
%> @date   12 May 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> Laplacian problem (quadrature-free implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyLaplacianQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for Laplacian problem (quadrature-free implementation).
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

function [Matrices] = VolumeMatricesAssemblyLaplacianQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        mu = Data.mu{id}(0,0);
        
        % Mass matrices assembly
        Matrices.Mprj_loc(1:nbases^2,1) = Integral.phiphiC;
    
        % Stiffness matrix assembly
        Matrices.A_loc(1:nbases^2,1) = mu*(Integral.gradxgradxC + Integral.gradygradyC);
 
end
                    