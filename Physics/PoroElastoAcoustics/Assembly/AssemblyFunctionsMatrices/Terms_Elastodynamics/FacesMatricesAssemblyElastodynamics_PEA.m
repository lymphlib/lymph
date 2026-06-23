%> @file   FacesMatricesAssemblyElastodynamics_PEA.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   5 June 2026
%> @brief Assembly of the matrices associated with face's integrals for 
%> elastodynamics problem  (preallocation and final assembly           
%> are associated with the functions in Physics/Elastodynamics).
%>
%==========================================================================
%> @section classFacesMatricesAssemblyElastodynamics_PEA Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for elastodynamics problem (preallocation and final assembly 
%> are associated with the functions in Physics/Elastodynamics).
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

function [Matrices] = FacesMatricesAssemblyElastodynamics_PEA(Data, femregion, Matrices, face, AssembInfo)
   
    ie = face.ie;
    iedg = face.iedg;
    bas_ie = 1:femregion.nbases(ie);

    % Evaluation of physical parameters
    vs    = Data.vs_el{femregion.id_phys(ie)}(face.xq,face.yq);
    vp    = Data.vp_el{femregion.id_phys(ie)}(face.xq,face.yq);
    mu    = Data.mu_el{femregion.id_phys(ie)}(face.xq,face.yq);
    lam   = Data.lam_el{femregion.id_phys(ie)}(face.xq,face.yq);
    
    % Internal face: Elastic neighbor
    if face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'E'
        mu_n   = Data.mu_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        lam_n  = Data.lam_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);

    % Internal face: Poroelastic neighbor
    elseif  face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'P'
        mu_n   = Data.mu{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        lam_n  = Data.lam{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);

    % Boundary face
    else
        mu_n   = Data.mu_el{femregion.id_phys(ie)}(face.xq,face.yq);
        lam_n  = Data.lam_el{femregion.id_phys(ie)}(face.xq,face.yq);
        
        % Absorbing boundary face (Parameters from I. Mazzieri PhD Thesis)
        if face.neigh_ie == -3
            c11 = (-mu ./ vs * face.ny * face.ny - (lam + 2*mu) ./ vp * face.nx * face.nx);
            c12 = ( mu ./ vs * face.nx * face.ny - (lam + 2*mu) ./ vp * face.nx * face.ny);
            c21 = ( mu ./ vs * face.nx * face.ny - (lam + 2*mu) ./ vp * face.nx * face.ny);
            c22 = (-mu ./ vs * face.nx * face.nx - (lam + 2*mu) ./ vp * face.ny * face.ny);
            c3_11 =  (mu .* (2*vs-vp) ./ vs + (lam .* vs + 2*mu .* (vs-vp)) ./ vp ) *face.nx*face.ny^2;
            c4_11 = -(mu .* (2*vs-vp) ./ vs + (lam .* vs + 2*mu .* (vs-vp)) ./ vp ) *face.ny*face.nx^2;
            c3_12 =  (mu .* (2*vs-vp) ./ vs * face.ny^3    - (lam .* vs + 2*mu .* (vs-vp)) ./ vp *face.ny*face.nx^2 );
            c4_12 = -(mu .* (2*vs-vp) ./ vs * face.nx*face.ny^2 + (lam .* vs + 2*mu .* (vs-vp)) ./ vp *face.nx^3 );
            c3_21 = -(mu .* (2*vs-vp) ./ vs * face.ny*face.nx^2 + (lam .* vs + 2*mu .* (vs-vp)) ./ vp *face.ny^3 );
            c4_21 =  (mu .* (2*vs-vp) ./ vs * face.nx^3    - (lam .* vs + 2*mu .* (vs-vp)) ./ vp *face.nx*face.ny^2 );
            c3_22 = -(mu .* (2*vs-vp) ./ vs + (lam .* vs + 2*mu .* (vs-vp)) ./ vp ) *face.nx*face.ny^2;
            c4_22 =  (mu .* (2*vs-vp) ./ vs + (lam .* vs + 2*mu .* (vs-vp)) ./ vp ) *face.ny*face.nx^2;
        end
    end

    lambda_ave = 2*lam .* lam_n ./ (lam + lam_n);
    mu_ave     = 2*mu .* mu_n ./ (mu + mu_n);
    harm_ave   = (lambda_ave + 2*mu_ave);            
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1
    
        Matrices.S1_P_loc(bas_ie,bas_ie) = Matrices.S1_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_P_loc(bas_ie,bas_ie) = Matrices.S4_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
    
        Matrices.IT1_loc(bas_ie,bas_ie) = Matrices.IT1_loc(bas_ie,bas_ie) + (face.ds .* (face.nx * (lambda_ave + 2*mu_ave) .* face.gradedgeqx + face.ny * mu_ave .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT2_loc(bas_ie,bas_ie) = Matrices.IT2_loc(bas_ie,bas_ie) + (face.ds .* (face.ny * lambda_ave .* face.gradedgeqx + face.nx * mu_ave .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT3_loc(bas_ie,bas_ie) = Matrices.IT3_loc(bas_ie,bas_ie) + (face.ds .* (face.ny * mu_ave .* face.gradedgeqx + face.nx * lambda_ave .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT4_loc(bas_ie,bas_ie) = Matrices.IT4_loc(bas_ie,bas_ie) + (face.ds .* (face.nx * mu_ave .* face.gradedgeqx + face.ny * (lambda_ave + 2*mu_ave) .* face.gradedgeqy ))' * face.phiedgeq;
    
    %% Absorbing boundary faces
    elseif face.neigh_ie == -3
    
        Matrices.ABC_S1_loc(bas_ie,bas_ie) = Matrices.ABC_S1_loc(bas_ie,bas_ie) + (face.ds .* c11 .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_S2_loc(bas_ie,bas_ie) = Matrices.ABC_S2_loc(bas_ie,bas_ie) + (face.ds .* c12 .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_S3_loc(bas_ie,bas_ie) = Matrices.ABC_S3_loc(bas_ie,bas_ie) + (face.ds .* c21 .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_S4_loc(bas_ie,bas_ie) = Matrices.ABC_S4_loc(bas_ie,bas_ie) + (face.ds .* c22 .* face.phiedgeq)' * face.phiedgeq;
    
        Matrices.ABC_R1_loc(bas_ie,bas_ie) = Matrices.ABC_R1_loc(bas_ie,bas_ie) + (face.ds .* (c3_11 .* face.gradedgeqx + c4_11 .* face.gradedgeqy))' * face.phiedgeq;
        Matrices.ABC_R2_loc(bas_ie,bas_ie) = Matrices.ABC_R2_loc(bas_ie,bas_ie) + (face.ds .* (c3_12 .* face.gradedgeqx + c4_12 .* face.gradedgeqy))' * face.phiedgeq;
        Matrices.ABC_R3_loc(bas_ie,bas_ie) = Matrices.ABC_R3_loc(bas_ie,bas_ie) + (face.ds .* (c3_21 .* face.gradedgeqx + c4_21 .* face.gradedgeqy))' * face.phiedgeq;
        Matrices.ABC_R4_loc(bas_ie,bas_ie) = Matrices.ABC_R4_loc(bas_ie,bas_ie) + (face.ds .* (c3_22 .* face.gradedgeqx + c4_22 .* face.gradedgeqy))' * face.phiedgeq;
    
    %% Internal faces: Elastic neighbor
    elseif face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'E'
    
        %% Element itself
        Matrices.S1_P_loc(bas_ie,bas_ie) = Matrices.S1_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_P_loc(bas_ie,bas_ie) = Matrices.S4_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
    
        Matrices.IT1_loc(bas_ie,bas_ie) = Matrices.IT1_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (face.nx * (lambda_ave + 2*mu_ave) .* face.gradedgeqx + face.ny * mu_ave .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT2_loc(bas_ie,bas_ie) = Matrices.IT2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (face.ny * lambda_ave .* face.gradedgeqx + face.nx * mu_ave .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT3_loc(bas_ie,bas_ie) = Matrices.IT3_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (face.ny * mu_ave .* face.gradedgeqx + face.nx * lambda_ave .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT4_loc(bas_ie,bas_ie) = Matrices.IT4_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (face.nx * mu_ave .* face.gradedgeqx + face.ny * (lambda_ave + 2*mu_ave) .* face.gradedgeqy ))' * face.phiedgeq;
    
        %% Neighboring element
        bas_neigh_ie = 1:femregion.nbases(face.neigh_ie);
        
        Matrices.S1_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S1_P_loc(face.neigh_idx,bas_neigh_ie) - (face.ds .* (harm_ave .* face.penalty_geom.max(iedg)) .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S4_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S4_P_loc(face.neigh_idx,bas_neigh_ie) - (face.ds .* (harm_ave .* face.penalty_geom.max(iedg)) .* face.phiedgeq)' * face.phiedgeqneigh;
    
        Matrices.IT1_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT1_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (face.nx * (lambda_ave + 2*mu_ave).* face.gradedgeqx + face.ny * mu_ave .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.IT2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (face.ny * lambda_ave .* face.gradedgeqx + face.nx * mu_ave .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.IT3_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT3_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (face.ny * mu_ave .* face.gradedgeqx + face.nx * lambda_ave .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.IT4_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT4_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (face.nx * mu_ave .* face.gradedgeqx + face.ny * (lambda_ave + 2*mu_ave) .* face.gradedgeqy ))' * face.phiedgeqneigh;
      
    %% Internal faces: Poroelastic neighbor
    elseif face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'P'
    
        %% Element itself
        Matrices.S1_P_loc(bas_ie,bas_ie) = Matrices.S1_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_P_loc(bas_ie,bas_ie) = Matrices.S4_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;

        Matrices.IT1_loc(bas_ie,bas_ie) = Matrices.IT1_loc(bas_ie,bas_ie) + (face.ds .* (face.nx * (lambda_ave + 2*mu_ave) .* face.phiedgeq))' *face.gradedgeqx + (face.ds .* (face.ny * mu_ave) .* face.phiedgeq)' * face.gradedgeqy;
        Matrices.IT2_loc(bas_ie,bas_ie) = Matrices.IT2_loc(bas_ie,bas_ie) + (face.ds .* (face.ny * lambda_ave .* face.phiedgeq))' * face.gradedgeqx + (face.ds .* (face.nx * mu_ave .* face.phiedgeq))' *face.gradedgeqy;
        Matrices.IT3_loc(bas_ie,bas_ie) = Matrices.IT3_loc(bas_ie,bas_ie) + (face.ds .* (face.ny * mu_ave .* face.phiedgeq))' *face.gradedgeqx + (face.ds .* (face.nx * lambda_ave .* face.phiedgeq))' * face.gradedgeqy;
        Matrices.IT4_loc(bas_ie,bas_ie) = Matrices.IT4_loc(bas_ie,bas_ie) + (face.ds .* (face.nx * mu_ave .*  face.phiedgeq))' * face.gradedgeqx + (face.ds .* (face.ny * (lambda_ave + 2*mu_ave) .* face.phiedgeq))' * face.gradedgeqy;
    
    end
     
end
