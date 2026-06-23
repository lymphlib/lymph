%> @file   FacesMatricesAssemblyHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   5 February 2026
%> @brief Assembly of the matrices associated with face's integrals for
%> heat equation.
%>
%==========================================================================
%> @section classFacesMatricesAssemblyHeat Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for heat equation.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Matrices struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = FacesMatricesAssemblyHeat(Data, femregion, Matrices, face, AssembInfo)
   
    ie = face.ie;
    iedg = face.iedg;

    % Evaluation of physical parameters
    mu = Data.mu{femregion.id(ie)}(face.xq,face.yq);
      
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        Matrices.IA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) = Matrices.IA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) + (face.ds .* mu .* (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy))' * face.phiedgeq;
        Matrices.SA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) = Matrices.SA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) + face.penalty_geom.max(iedg) * (face.ds .* mu .* face.phiedgeq)' * face.phiedgeq;
    
    %% Internal faces
    elseif face.neigh_ie > 0

        %% Element itself
        Matrices.IA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) = Matrices.IA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) + 0.5 * (face.ds .* mu .* (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy))' * face.phiedgeq;
        Matrices.SA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) = Matrices.SA_loc(1:femregion.nbases(ie),1:femregion.nbases(ie)) + face.penalty_geom.max(iedg) * (face.ds .* mu .* face.phiedgeq)' * face.phiedgeq;

        %% Neighboring element
        Matrices.IA_loc(face.neigh_idx,1:femregion.nbases(face.neigh_ie)) = Matrices.IA_loc(face.neigh_idx,1:femregion.nbases(face.neigh_ie)) - 0.5 * (face.ds .* mu .* ( face.nx * face.gradedgeqx +  face.ny * face.gradedgeqy))' * face.phiedgeqneigh;
        Matrices.SA_loc(face.neigh_idx,1:femregion.nbases(face.neigh_ie)) = Matrices.SA_loc(face.neigh_idx,1:femregion.nbases(face.neigh_ie)) - face.penalty_geom.max(iedg) * (face.ds .* mu .* face.phiedgeq)' * face.phiedgeqneigh;

    end

end