%> @file   VolumeMatricesAssemblyAcousticsQF.m
%> @author Mattia Corti
%> @date   5 June 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> acoustics problem (quadrature-free implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyAcousticsQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for acoustics problem (quadrature-free implementation).
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

function [Matrices] = VolumeMatricesAssemblyAcousticsQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        rho_a   = Data.rho_a{id}(0,0);
        c       = Data.c{id}(0,0);
                    
        % Matrices assembly
        Matrices.W_loc(1:nbases^2) = rho_a*(Integral.gradxgradxC+Integral.gradygradyC);

        Matrices.M_P_A_loc(1:nbases^2) = (c.^(-2).* rho_a)*Integral.phiphiC;
        Matrices.MPrjA_loc(1:nbases^2) = Integral.phiphiC;
 
end
                    
