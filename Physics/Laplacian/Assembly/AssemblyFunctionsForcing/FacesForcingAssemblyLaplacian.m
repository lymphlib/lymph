%> @file   FacesForcingAssemblyLaplacian.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   12 May 2026
%> @brief Assembly of the forcing vectors associated with face's integrals for
%> Laplacian problem.
%>
%==========================================================================
%> @section classFacesForcingAssemblyLaplacian Class description
%==========================================================================
%> @brief            Assembly of the forcing vectors associated with face's
%> integrals for Laplacian problem.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Forcing    Forcing struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = FacesForcingAssemblyLaplacian(Data, femregion, Forcing, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;

    % Evaluation of physical parameters
    mu = Data.mu{femregion.id(ie)}(face.xq,face.yq);    
    gD = Data.DirBC{femregion.id(ie)}(face.xq,face.yq);
     
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        Forcing.F_D_loc(1:femregion.nbases(ie),1) = Forcing.F_D_loc(1:femregion.nbases(ie),1) - (face.ds .* mu .* ( face.nx * face.gradedgeqx + face.ny * face.gradedgeqy))' * gD;
        Forcing.F_D_loc(1:femregion.nbases(ie),1) = Forcing.F_D_loc(1:femregion.nbases(ie),1) + face.penalty_geom.max(iedg) * (face.ds .* mu .* face.phiedgeq)' * gD;        
        
    end

end
