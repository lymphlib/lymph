%> @file   IPVolumeMatricesAssemblyFHNST.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief Assembly of the matrices associated with volume integrals of bilinear forms for FHN
%> problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classIPVolumeMatricesAssemblyFHNST Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals of bilinear forms for FHN problem (subtriangulation implementation).
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

function [Matrices] = IPVolumeMatricesAssemblyFHNST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        D_ext  = Data.D_ext{id}(elem.xq,elem.yq,0);
        Chi    = Data.Chi{id}(elem.xq,elem.yq,0);
        Cm     = Data.Cm{id}(elem.xq,elem.yq,0);

        if ~Data.isotropy
            D_axn  = Data.D_axn{id}(elem.xq,elem.xq,0);

            % Axonal directions
            D_axn_xx = AssembInfo.Tens.D_xx(ie);
            D_axn_xy = AssembInfo.Tens.D_xy(ie);
            D_axn_yx = AssembInfo.Tens.D_yx(ie);
            D_axn_yy = AssembInfo.Tens.D_yy(ie);
        end

        % Mass matrices assembly
        Matrices.M_prj_loc(1:nbases,1:nbases)   = (elem.dx .* elem.phiq)' * elem.phiq;
        Matrices.M_u_loc(1:nbases,1:nbases)     = (elem.dx .* Cm .* Chi .* elem.phiq)' * elem.phiq;
        Matrices.M_w_loc(1:nbases,1:nbases)     = (elem.dx .*  elem.phiq)' * elem.phiq;
              
        % Stiffness matrix assembly
        Matrices.A_loc(1:nbases,1:nbases) = (elem.dx .* (D_ext .* elem.gradqx))' * elem.gradqx ...
                                              + (elem.dx .* (D_ext .* elem.gradqy))' * elem.gradqy;

        if ~Data.isotropy
            Matrices.A_loc(1:nbases,1:nbases) = Matrices.A_loc(1:nbases,1:nbases) ...
                                                  + (elem.dx .* (D_axn .* D_axn_xx .* elem.gradqx))' * elem.gradqx ...
                                                  + (elem.dx .* (D_axn .* D_axn_yx .* elem.gradqy))' * elem.gradqx ...
                                                  + (elem.dx .* (D_axn .* D_axn_xy .* elem.gradqx))' * elem.gradqy ...
                                                  + (elem.dx .* (D_axn .* D_axn_yy .* elem.gradqy))' * elem.gradqy;
        end

end
                    