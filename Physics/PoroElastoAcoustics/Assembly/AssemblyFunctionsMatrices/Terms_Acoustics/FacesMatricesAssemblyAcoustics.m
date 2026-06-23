%> @file   FacesMatricesAssemblyAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Assembly of the matrices associated with face's integrals for
%> elastodynamics problem.
%>
%==========================================================================
%> @section classFacesMatricesAssemblyAcoustics Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for acoustics problem.
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

function [Matrices] = FacesMatricesAssemblyAcoustics(Data, femregion, Matrices, face, AssembInfo)
   
    ie = face.ie;
    iedg = face.iedg;
    bas_ie = 1:femregion.nbases(ie);

    % Evaluation of physical parameters
    rho_a   = Data.rho_a{femregion.id_phys(ie)}(face.xq,face.yq);
    
    % Acoustic neighbor
    if face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'A'
        rho_a_n   = Data.rho_a{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        rho_a_ave = 2*rho_a.*rho_a_n./(rho_a + rho_a_n);
    
    % Boundary edge
    elseif face.neigh_ie < 0
        rho_a_n   = Data.rho_a{femregion.id_phys(ie)}(face.xq,face.yq);
        rho_a_ave = 2*rho_a.*rho_a_n./(rho_a + rho_a_n);
    
    end
    
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1
    
        Matrices.S_A_loc(bas_ie,bas_ie)  = Matrices.S_A_loc(bas_ie,bas_ie) + face.penalty_geom.harm(iedg) * (face.ds .* rho_a_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.IT_A_loc(bas_ie,bas_ie) = Matrices.IT_A_loc(bas_ie,bas_ie) + (face.ds .* (face.nx * rho_a_ave .* face.gradedgeqx + face.ny * rho_a_ave .* face.gradedgeqy))' * face.phiedgeq;

    %% Internal faces
    elseif face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'A'
    
        %% Element itself
        Matrices.S_A_loc(bas_ie,bas_ie)  = Matrices.S_A_loc(bas_ie,bas_ie) + face.penalty_geom.harm(iedg) * (face.ds .* rho_a_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.IT_A_loc(bas_ie,bas_ie) = Matrices.IT_A_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (face.nx * rho_a_ave .* face.gradedgeqx + face.ny * rho_a_ave .* face.gradedgeqy))' * face.phiedgeq;

        %% Neighboring element
        bas_neigh_ie = 1:femregion.nbases(face.neigh_ie);
        
        Matrices.S_A_loc(face.neigh_idx,bas_neigh_ie)  = Matrices.S_A_loc(face.neigh_idx,bas_neigh_ie) - (face.ds .* (rho_a_ave .* face.penalty_geom.harm(iedg)) .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.IT_A_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT_A_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (face.nx * rho_a_ave .* face.gradedgeqx + face.ny * rho_a_ave .* face.gradedgeqy))' * face.phiedgeqneigh;

    end
     
end