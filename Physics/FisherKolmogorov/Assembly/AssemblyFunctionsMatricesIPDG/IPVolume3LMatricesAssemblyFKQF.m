%> @file   IPVolume3LMatricesAssemblyFKQF.m
%> @author Mattia Corti
%> @date   16 September 2025
%> @brief Assembly of the matrices associated with volume integrals of trilinear forms for FK
%> problem (quadrature-free implementation).
%>
%==========================================================================
%> @section classIPVolume3LMatricesAssemblyFKQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals of trilinear forms for FK problem (quadrature-free implementation).
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

function [Matrices] = IPVolume3LMatricesAssemblyFKQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        alpha  = Data.alpha{id}(0,0,0);

        % 3L mass matrix assembly
        Matrices.M_NL_loc(1:nbases^2,1:nbases) = reshape(alpha*Integral.phiphiphiC,[nbases^2, nbases]);
      
end
                    
