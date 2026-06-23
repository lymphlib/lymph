%> @file   FacesForcingAssemblyFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief Assembly of the Dirichlet BCs face integrals for FHN problem.
%>
%==========================================================================
%> @section classFacesForcingAssemblyFHN Class description
%==========================================================================
%> @brief            Assembly of the Dirichlet BCs face integrals for FHN problem.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Forcing    Forcing struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the assembly procedure
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = FacesForcingAssemblyFHN(Data, femregion, Forcing, face, AssembInfo)
   
        ie = face.ie;
        iedg = face.iedg;

        % Evaluation of physical parameters
        D_ext  = Data.D_ext{femregion.id(ie)}(face.xq,face.yq,AssembInfo.t);
        gD     = Data.DirBC{femregion.id(ie)}(face.xq,face.yq,AssembInfo.t);
     
        if ~Data.isotropy
            D_axn  = Data.D_axn{femregion.id(ie)}(face.xq,face.yq,AssembInfo.t);

            % Axonal directions
            D_axn_xx = AssembInfo.Tens.D_xx(ie);
            D_axn_xy = AssembInfo.Tens.D_xy(ie);
            D_axn_yx = AssembInfo.Tens.D_yx(ie);
            D_axn_yy = AssembInfo.Tens.D_yy(ie);
        end
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        if ~Data.isotropy
            penalty = (D_ext+D_axn) .* face.penalty_geom.harm(iedg);
        else
            penalty = D_ext .* face.penalty_geom.harm(iedg);
        end

        % Extracellular component
        Forcing.F_D_loc(1:femregion.nbases(ie),1) = Forcing.F_D_loc(1:femregion.nbases(ie),1) - (face.ds .* D_ext .* (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy))' * gD;
        Forcing.F_D_loc(1:femregion.nbases(ie),1) = Forcing.F_D_loc(1:femregion.nbases(ie),1) + (face.ds .* penalty .* face.phiedgeq)' * gD;        
                 
        % Axonal Component
        if ~Data.isotropy
            Forcing.F_D_loc(1:femregion.nbases(ie),1) = Forcing.F_D_loc(1:femregion.nbases(ie),1) - (face.ds .* D_axn * (D_axn_xx * face.nx + D_axn_yx * face.ny) .* face.gradedgeqx)' * gD;
            Forcing.F_D_loc(1:femregion.nbases(ie),1) = Forcing.F_D_loc(1:femregion.nbases(ie),1) - (face.ds .* D_axn * (D_axn_xy * face.nx + D_axn_yy * face.ny) .* face.gradedgeqy)' * gD;
        end
    
    end

end
