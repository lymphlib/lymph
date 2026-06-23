%> @file   IPVolumeMatricesAssemblyFKQF.m
%> @author Mattia Corti
%> @date   14 April 2026
%> @brief Assembly of the matrices associated with volume integrals of bilinear forms for FK
%> problem (quadrature-free implementation).
%>
%==========================================================================
%> @section classIPVolumeMatricesAssemblyFKQF Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals of bilinear forms for FK problem (quadrature-free implementation).
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

function [Matrices] = IPVolumeMatricesAssemblyFKQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        D_ext  = Data.D_ext{id}(0,0,0);
        alpha  = Data.alpha{id}(0,0,0);

        if ~Data.isotropy
            D_axn  = Data.D_axn{id}(0,0,0);
            
            % Axonal directions
            D_axn_xx = AssembInfo.Tens.D_xx(ie);
            D_axn_xy = AssembInfo.Tens.D_xy(ie);
            D_axn_yx = AssembInfo.Tens.D_yx(ie);
            D_axn_yy = AssembInfo.Tens.D_yy(ie);
        end

        % Mass matrices assembly
        Matrices.M_prj_loc(1:nbases^2,1)  = Integral.phiphiC;
        Matrices.M_loc(1:nbases^2,1)      = alpha*Integral.phiphiC;
       
        % Stiffness matrix assembly
        Matrices.A_loc(1:nbases^2,1) = D_ext*(Integral.gradxgradxC + Integral.gradygradyC);

        if ~Data.isotropy
            Matrices.A_loc(1:nbases^2,1) = Matrices.A_loc(1:nbases^2,ie) ...
                                            + D_axn*(D_axn_xx*Integral.gradxgradxC ...
                                            + D_axn_xy*Integral.gradxgradyC ...
                                            + D_axn_yx*Integral.gradygradxC ...
                                            + D_axn_yy*Integral.gradygradyC);
        end
   
end
                    