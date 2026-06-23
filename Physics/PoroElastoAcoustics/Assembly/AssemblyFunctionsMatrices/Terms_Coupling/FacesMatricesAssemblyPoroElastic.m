%> @file   FacesMatricesAssemblyPoroElastic.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   5 June 2026
%> @brief Assembly of the matrices associated with face's integrals for
%> poroelasticity-elastodynamics coupling.
%>
%==========================================================================
%> @section classFacesMatricesAssemblyPoroElastic Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for poroelasticity-elastodynamics coupling.
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

function [Matrices] = FacesMatricesAssemblyPoroElastic(Data, femregion, Matrices, face, AssembInfo)
   
    ie = face.ie;
    iedg = face.iedg;
    bas_ie = 1:femregion.nbases(ie);

    % Evaluation of physical parameters
    mu     = Data.mu{femregion.id_phys(ie)}(face.xq,face.yq);
    lam    = Data.lam{femregion.id_phys(ie)}(face.xq,face.yq);
    mu_n   = Data.mu_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
    lam_n  = Data.lam_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);

    % Auxiliary quantities (cf. physical parameters) for vector assembling
    lam_ave  = 2 * lam .* lam_n ./ (lam + lam_n);
    mu_ave   = 2 * mu .* mu_n ./ (mu + mu_n);
    harm_ave = lam_ave + 2*mu_ave;
    
    bas_neigh_ie = 1:femregion.nbases(face.neigh_ie);
        
    Matrices.C1_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C1_P_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* (lam_ave + 2*mu_ave) .* face.nx .* face.phiedgeq)' * face.gradedgeqxneigh ...
                                                    - ( face.ds .* mu_ave .* face.ny .* face.phiedgeq)' * face.gradedgeqyneigh ...
                                                    - face.penalty_geom.max(iedg) * ( face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeqneigh;

    Matrices.C2_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C2_P_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* lam_ave .* face.nx .* face.phiedgeq)' * face.gradedgeqyneigh ...
                                                    - ( face.ds .* mu_ave .* face.ny .* face.phiedgeq)' * face.gradedgeqxneigh;

    Matrices.C3_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C3_P_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* mu_ave .* face.nx .* face.phiedgeq)' * face.gradedgeqyneigh ...
                                                    - ( face.ds .* lam_ave .* face.ny .* face.phiedgeq)' * face.gradedgeqxneigh;

    Matrices.C4_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C4_P_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* mu_ave .* face.nx .* face.phiedgeq)' * face.gradedgeqxneigh ...
                                                    - ( face.ds .* (lam_ave + 2*mu_ave) .* face.ny .* face.phiedgeq)' * face.gradedgeqyneigh ...
                                                    - face.penalty_geom.max(iedg) * ( face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeqneigh;

    Matrices.C1_E_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C1_E_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* (lam_ave + 2*mu_ave) .* face.nx .* face.phiedgeq)' * face.gradedgeqxneigh ...
                                                    - ( face.ds .* mu_ave .* face.ny .* face.phiedgeq)' *  face.gradedgeqyneigh ...
                                                    - face.penalty_geom.max(iedg) * ( face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeqneigh;

    Matrices.C2_E_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C2_E_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* mu_ave .* face.nx .* face.phiedgeq)' * face.gradedgeqyneigh ...
                                                    - ( face.ds .* lam_ave .* face.ny .* face.phiedgeq)' * face.gradedgeqxneigh;

    Matrices.C3_E_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C3_E_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* lam_ave .* face.nx .* face.phiedgeq)' * face.gradedgeqyneigh ...
                                                    - ( face.ds .* mu_ave .* face.ny .* face.phiedgeq)' * face.gradedgeqxneigh;

    Matrices.C4_E_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C4_E_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    - ( face.ds .* mu_ave .* face.nx .* face.phiedgeq)' * face.gradedgeqxneigh ...
                                                    - ( face.ds .* (lam_ave + 2*mu_ave) .* face.ny .* face.phiedgeq)' *  face.gradedgeqyneigh ...
                                                    - face.penalty_geom.max(iedg) * ( face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeqneigh;
                                                 
end
