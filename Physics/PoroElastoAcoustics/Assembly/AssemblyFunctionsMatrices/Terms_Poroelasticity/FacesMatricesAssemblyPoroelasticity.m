%> @file   FacesMatricesAssemblyPoroelasticity.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   5 June 2026
%> @brief Assembly of the matrices associated with face's integrals for
%> poroelasticity problem.
%>
%==========================================================================
%> @section classFacesMatricesAssemblyPoroelasticity Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for poroelasticity problem.
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

function [Matrices] = FacesMatricesAssemblyPoroelasticity(Data, femregion, Matrices, face, AssembInfo)
   
    ie = face.ie;
    iedg = face.iedg;
    bas_ie = 1:femregion.nbases(ie);

    % Evaluation of physical parameters
    mu        = Data.mu{femregion.id_phys(ie)}(face.xq,face.yq);
    lam       = Data.lam{femregion.id_phys(ie)}(face.xq,face.yq);
    beta      = Data.beta{femregion.id_phys(ie)}(face.xq,face.yq);
    m         = Data.m{femregion.id_phys(ie)}(face.xq,face.yq);

    % Boundary edge
    if face.neigh_ie < 0
        mu_n      = Data.mu{femregion.id_phys(ie)}(face.xq,face.yq);
        lam_n     = Data.lam{femregion.id_phys(ie)}(face.xq,face.yq);
        m_n       = Data.m{femregion.id_phys(ie)}(face.xq,face.yq);
        beta_n    = Data.beta{femregion.id_phys(ie)}(face.xq,face.yq);
        vp_poroI  = Data.vp_poroI{femregion.id_phys(ie)}(face.xq,face.yq);
        vp_poroII = Data.vp_poroII{femregion.id_phys(ie)}(face.xq,face.yq);
        vs_poro   = Data.vs_poro{femregion.id_phys(ie)}(face.xq,face.yq);
        phi_por   = Data.phi_por{femregion.id_phys(ie)}(face.xq,face.yq);
        a_coef    = Data.a_coef{femregion.id_phys(ie)}(face.xq,face.yq);

    % Elastic neighbor
    elseif femregion.label(face.neigh_ie) == 'E'
        mu_n   = Data.mu_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        lam_n  = Data.lam_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        m_n    = Data.m{femregion.id_phys(ie)}(face.xq,face.yq);
        beta_n = Data.beta{femregion.id_phys(ie)}(face.xq,face.yq);
    
    % Poroelastic neighbor
    elseif femregion.label(face.neigh_ie) == 'P'
        mu_n   = Data.mu{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        lam_n  = Data.lam{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        m_n    = Data.m{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
        beta_n = Data.beta{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);

    % Acoustic neighbor
    elseif femregion.label(face.neigh_ie) == 'A'
        mu_n      = Data.mu_el{femregion.id_phys(ie)}(face.xq,face.yq);
        lam_n     = Data.lam_el{femregion.id_phys(ie)}(face.xq,face.yq);
        m_n       = Data.m{femregion.id_phys(ie)}(face.xq,face.yq);
        beta_n    = Data.beta{femregion.id_phys(ie)}(face.xq,face.yq);

    end

    lambda_ave = 2*lam .* lam_n ./ (lam + lam_n);
    mu_ave     = 2*mu .* mu_n ./ (mu + mu_n);
    harm_ave   = (lambda_ave + 2*mu_ave);
    m_ave      = 2*m .* m_n ./ (m + m_n);
    beta_ave   = 2*beta .* beta_n ./ (beta + beta_n);
           
    %% Dirichlet boundary faces
    if face.neigh_ie == -1
    
        Matrices.S1_P_loc(bas_ie,bas_ie) = Matrices.S1_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_P_loc(bas_ie,bas_ie) = Matrices.S4_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_loc(bas_ie,bas_ie) = Matrices.S1_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_loc(bas_ie,bas_ie) = Matrices.S2_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_loc(bas_ie,bas_ie) = Matrices.S3_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_loc(bas_ie,bas_ie) = Matrices.S4_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_beta_loc(bas_ie,bas_ie) = Matrices.S1_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_beta_loc(bas_ie,bas_ie) = Matrices.S2_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_beta_loc(bas_ie,bas_ie) = Matrices.S3_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_beta_loc(bas_ie,bas_ie) = Matrices.S4_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_beta2_loc(bas_ie,bas_ie) = Matrices.S1_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_beta2_loc(bas_ie,bas_ie) = Matrices.S2_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_beta2_loc(bas_ie,bas_ie) = Matrices.S3_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_beta2_loc(bas_ie,bas_ie) = Matrices.S4_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.IT1_loc(bas_ie,bas_ie) = Matrices.IT1_loc(bas_ie,bas_ie) + (face.ds .* ((lambda_ave + 2*mu_ave) .* face.nx .* face.gradedgeqx + mu_ave .* face.ny .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT2_loc(bas_ie,bas_ie) = Matrices.IT2_loc(bas_ie,bas_ie) + (face.ds .* (lambda_ave .* face.ny .* face.gradedgeqx + mu_ave .* face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT3_loc(bas_ie,bas_ie) = Matrices.IT3_loc(bas_ie,bas_ie) + (face.ds .* (mu_ave .* face.ny .* face.gradedgeqx + lambda_ave .* face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT4_loc(bas_ie,bas_ie) = Matrices.IT4_loc(bas_ie,bas_ie) + (face.ds .* (mu_ave .* face.nx .* face.gradedgeqx + (lambda_ave + 2*mu_ave) .* face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        Matrices.BT1_loc(bas_ie,bas_ie) = Matrices.BT1_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .*face.nx .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT2_loc(bas_ie,bas_ie) = Matrices.BT2_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .*face.ny .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT3_loc(bas_ie,bas_ie) = Matrices.BT3_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .*face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.BT4_loc(bas_ie,bas_ie) = Matrices.BT4_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .*face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        Matrices.BT1_beta_loc(bas_ie,bas_ie) = Matrices.BT1_beta_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave .*face.nx .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT2_beta_loc(bas_ie,bas_ie) = Matrices.BT2_beta_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave .*face.ny .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT3_beta_loc(bas_ie,bas_ie) = Matrices.BT3_beta_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave .*face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.BT4_beta_loc(bas_ie,bas_ie) = Matrices.BT4_beta_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave .*face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        Matrices.BT1_beta2_loc(bas_ie,bas_ie) = Matrices.BT1_beta2_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave.^2 .*face.nx .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT2_beta2_loc(bas_ie,bas_ie) = Matrices.BT2_beta2_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave.^2 .*face.ny .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT3_beta2_loc(bas_ie,bas_ie) = Matrices.BT3_beta2_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave.^2 .*face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.BT4_beta2_loc(bas_ie,bas_ie) = Matrices.BT4_beta2_loc(bas_ie,bas_ie) + (face.ds .* (m_ave .* beta_ave.^2 .*face.ny .* face.gradedgeqy ))' * face.phiedgeq;

    %% Absorbing boundary faces
    elseif face.neigh_ie == -3
    
        Matrices.ABC_uu_1_loc(bas_ie,bas_ie) = Matrices.ABC_uu_1_loc(bas_ie,bas_ie) + (face.ds .* (rho.*vp_poroI*face.nx*face.nx + (rho - rho_f.*phi_por./a_coef).*vs_poro*(1-face.nx*face.nx)) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_uu_2_loc(bas_ie,bas_ie) = Matrices.ABC_uu_2_loc(bas_ie,bas_ie) + (face.ds .* (rho.*vp_poroI*face.nx*face.ny - (rho - rho_f.*phi_por./a_coef).*vs_poro*face.nx*face.ny)     .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_uu_3_loc(bas_ie,bas_ie) = Matrices.ABC_uu_3_loc(bas_ie,bas_ie) + (face.ds .* (rho.*vp_poroI*face.ny*face.nx - (rho - rho_f.*phi_por./a_coef).*vs_poro*face.ny*face.nx)     .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_uu_4_loc(bas_ie,bas_ie) = Matrices.ABC_uu_4_loc(bas_ie,bas_ie) + (face.ds .* (rho.*vp_poroI*face.ny*face.ny + (rho - rho_f.*phi_por./a_coef).*vs_poro*(1-face.ny*face.ny)) .* face.phiedgeq)' * face.phiedgeq;

        Matrices.ABC_uw_1_loc(bas_ie,bas_ie) = Matrices.ABC_uu_1_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroII*face.nx*face.nx) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_uw_2_loc(bas_ie,bas_ie) = Matrices.ABC_uu_2_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroII*face.nx*face.ny) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_uw_3_loc(bas_ie,bas_ie) = Matrices.ABC_uu_3_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroII*face.ny*face.nx) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_uw_4_loc(bas_ie,bas_ie) = Matrices.ABC_uu_4_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroII*face.ny*face.ny) .* face.phiedgeq)' * face.phiedgeq;

        Matrices.ABC_wu_1_loc(bas_ie,bas_ie) = Matrices.ABC_uu_1_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroI*face.nx*face.nx) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_wu_2_loc(bas_ie,bas_ie) = Matrices.ABC_uu_2_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroI*face.nx*face.ny) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_wu_3_loc(bas_ie,bas_ie) = Matrices.ABC_uu_3_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroI*face.ny*face.nx) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_wu_4_loc(bas_ie,bas_ie) = Matrices.ABC_uu_4_loc(bas_ie,bas_ie) + (face.ds .* (rho_f.*vp_poroI*face.ny*face.ny) .* face.phiedgeq)' * face.phiedgeq;

        Matrices.ABC_ww_1_loc(bas_ie,bas_ie) = Matrices.ABC_uu_1_loc(bas_ie,bas_ie) + (face.ds .* (rho_w.*vp_poroII*face.nx*face.nx) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_ww_2_loc(bas_ie,bas_ie) = Matrices.ABC_uu_2_loc(bas_ie,bas_ie) + (face.ds .* (rho_w.*vp_poroII*face.nx*face.ny) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_ww_3_loc(bas_ie,bas_ie) = Matrices.ABC_uu_3_loc(bas_ie,bas_ie) + (face.ds .* (rho_w.*vp_poroII*face.ny*face.nx) .* face.phiedgeq)' * face.phiedgeq;
        Matrices.ABC_ww_4_loc(bas_ie,bas_ie) = Matrices.ABC_uu_4_loc(bas_ie,bas_ie) + (face.ds .* (rho_w.*vp_poroII*face.ny*face.ny) .* face.phiedgeq)' * face.phiedgeq;

    %% Internal faces with poroelastic neighbor
    elseif face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'P'
    
        %% Element itself
        Matrices.S1_P_loc(bas_ie,bas_ie) = Matrices.S1_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_P_loc(bas_ie,bas_ie) = Matrices.S4_P_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_loc(bas_ie,bas_ie) = Matrices.S1_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_loc(bas_ie,bas_ie) = Matrices.S2_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_loc(bas_ie,bas_ie) = Matrices.S3_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_loc(bas_ie,bas_ie) = Matrices.S4_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_beta_loc(bas_ie,bas_ie) = Matrices.S1_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_beta_loc(bas_ie,bas_ie) = Matrices.S2_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_beta_loc(bas_ie,bas_ie) = Matrices.S3_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_beta_loc(bas_ie,bas_ie) = Matrices.S4_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_beta2_loc(bas_ie,bas_ie) = Matrices.S1_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_beta2_loc(bas_ie,bas_ie) = Matrices.S2_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_beta2_loc(bas_ie,bas_ie) = Matrices.S3_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_beta2_loc(bas_ie,bas_ie) = Matrices.S4_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.IT1_loc(bas_ie,bas_ie) = Matrices.IT1_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* ((lambda_ave + 2*mu_ave) .* face.nx .* face.gradedgeqx + mu_ave .* face.ny .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT2_loc(bas_ie,bas_ie) = Matrices.IT2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (lambda_ave .* face.ny .* face.gradedgeqx + mu_ave .* face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT3_loc(bas_ie,bas_ie) = Matrices.IT3_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (mu_ave .* face.ny .* face.gradedgeqx + lambda_ave .* face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.IT4_loc(bas_ie,bas_ie) = Matrices.IT4_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (mu_ave .* face.nx .* face.gradedgeqx + (lambda_ave + 2*mu_ave) .* face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        Matrices.BT1_loc(bas_ie,bas_ie) = Matrices.BT1_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .*face.nx .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT2_loc(bas_ie,bas_ie) = Matrices.BT2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .*face.ny .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT3_loc(bas_ie,bas_ie) = Matrices.BT3_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .*face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.BT4_loc(bas_ie,bas_ie) = Matrices.BT4_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .*face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        Matrices.BT1_beta_loc(bas_ie,bas_ie) = Matrices.BT1_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave .*face.nx .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT2_beta_loc(bas_ie,bas_ie) = Matrices.BT2_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave .*face.ny .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT3_beta_loc(bas_ie,bas_ie) = Matrices.BT3_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave .*face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.BT4_beta_loc(bas_ie,bas_ie) = Matrices.BT4_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave .*face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        Matrices.BT1_beta2_loc(bas_ie,bas_ie) = Matrices.BT1_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.nx .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT2_beta2_loc(bas_ie,bas_ie) = Matrices.BT2_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.ny .* face.gradedgeqx ))' * face.phiedgeq;
        Matrices.BT3_beta2_loc(bas_ie,bas_ie) = Matrices.BT3_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.nx .* face.gradedgeqy ))' * face.phiedgeq;
        Matrices.BT4_beta2_loc(bas_ie,bas_ie) = Matrices.BT4_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.ny .* face.gradedgeqy ))' * face.phiedgeq;

        %% Neighboring element
        bas_neigh_ie = 1:femregion.nbases(face.neigh_ie);
        
        Matrices.S1_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S1_P_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S4_P_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S4_P_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeqneigh;

        Matrices.S1_B_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S1_B_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S2_B_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S2_B_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S3_B_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S3_B_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S4_B_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S4_B_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;

        Matrices.S1_B_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S1_B_beta_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S2_B_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S2_B_beta_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S3_B_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S3_B_beta_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S4_B_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S4_B_beta_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;

        Matrices.S1_B_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S1_B_beta2_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S2_B_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S2_B_beta2_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S3_B_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S3_B_beta2_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;
        Matrices.S4_B_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.S4_B_beta2_loc(face.neigh_idx,bas_neigh_ie) - face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;

        Matrices.IT1_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT1_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* ((lambda_ave + 2*mu_ave) .* face.nx .* face.gradedgeqx + mu_ave .* face.ny .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.IT2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (lambda_ave .* face.ny .* face.gradedgeqx + mu_ave .* face.nx .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.IT3_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT3_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (mu_ave .* face.ny .* face.gradedgeqx + lambda_ave .* face.nx .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.IT4_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IT4_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (mu_ave .* face.nx .* face.gradedgeqx + (lambda_ave + 2*mu_ave) .* face.ny .* face.gradedgeqy ))' * face.phiedgeqneigh;

        Matrices.BT1_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT1_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .*face.nx .* face.gradedgeqx ))' * face.phiedgeqneigh;
        Matrices.BT2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .*face.ny .* face.gradedgeqx ))' * face.phiedgeqneigh;
        Matrices.BT3_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT3_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .*face.nx .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.BT4_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT4_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .*face.ny .* face.gradedgeqy ))' * face.phiedgeqneigh;

        Matrices.BT1_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT1_beta_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave .*face.nx .* face.gradedgeqx ))' * face.phiedgeqneigh;
        Matrices.BT2_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT2_beta_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave .*face.ny .* face.gradedgeqx ))' * face.phiedgeqneigh;
        Matrices.BT3_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT3_beta_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave .*face.nx .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.BT4_beta_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT4_beta_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave .*face.ny .* face.gradedgeqy ))' * face.phiedgeqneigh;

        Matrices.BT1_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT1_beta2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.nx .* face.gradedgeqx ))' * face.phiedgeqneigh;
        Matrices.BT2_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT2_beta2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.ny .* face.gradedgeqx ))' * face.phiedgeqneigh;
        Matrices.BT3_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT3_beta2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.nx .* face.gradedgeqy ))' * face.phiedgeqneigh;
        Matrices.BT4_beta2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.BT4_beta2_loc(face.neigh_idx,bas_neigh_ie) - 0.5 * (face.ds .* (m_ave .* beta_ave.^2 .*face.ny .* face.gradedgeqy ))' * face.phiedgeqneigh;

    %% Internal faces with elastic neighbor
    elseif face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'E'
    
        %% Element itself
        Matrices.S1_P_EP_loc(bas_ie,bas_ie) = Matrices.S1_P_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_P_EP_loc(bas_ie,bas_ie) = Matrices.S4_P_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* harm_ave .* face.phiedgeq)' * face.phiedgeq;

        Matrices.BT1_beta2_EP_loc(bas_ie,bas_ie) = Matrices.BT1_beta2_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.phiedgeq)' * face.gradedgeqx;
        Matrices.BT2_beta2_EP_loc(bas_ie,bas_ie) = Matrices.BT2_beta2_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.phiedgeq)' * face.gradedgeqy;
        Matrices.BT3_beta2_EP_loc(bas_ie,bas_ie) = Matrices.BT3_beta2_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.phiedgeq)' * face.gradedgeqx;
        Matrices.BT4_beta2_EP_loc(bas_ie,bas_ie) = Matrices.BT4_beta2_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.phiedgeq)' * face.gradedgeqy;

        Matrices.S1_B_beta2_EP_loc(bas_ie,bas_ie) = Matrices.S1_B_beta2_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_beta2_EP_loc(bas_ie,bas_ie) = Matrices.S2_B_beta2_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_beta2_EP_loc(bas_ie,bas_ie) = Matrices.S3_B_beta2_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_beta2_EP_loc(bas_ie,bas_ie) = Matrices.S4_B_beta2_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.S1_B_beta_EP_loc(bas_ie,bas_ie) = Matrices.S1_B_beta_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_beta_EP_loc(bas_ie,bas_ie) = Matrices.S2_B_beta_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_beta_EP_loc(bas_ie,bas_ie) = Matrices.S3_B_beta_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_beta_EP_loc(bas_ie,bas_ie) = Matrices.S4_B_beta_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.BT1_beta_EP_loc(bas_ie,bas_ie) = Matrices.BT1_beta_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave .* face.nx .* face.gradedgeqx)' *  face.phiedgeq;
        Matrices.BT2_beta_EP_loc(bas_ie,bas_ie) = Matrices.BT2_beta_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave .* face.ny .* face.gradedgeqx)' *  face.phiedgeq;
        Matrices.BT3_beta_EP_loc(bas_ie,bas_ie) = Matrices.BT3_beta_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave .* face.nx .* face.gradedgeqy)' *  face.phiedgeq;
        Matrices.BT4_beta_EP_loc(bas_ie,bas_ie) = Matrices.BT4_beta_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* beta_ave .* face.ny .* face.gradedgeqy)' *  face.phiedgeq;

        Matrices.S1_B_EP_loc(bas_ie,bas_ie) = Matrices.S1_B_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S2_B_EP_loc(bas_ie,bas_ie) = Matrices.S2_B_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S3_B_EP_loc(bas_ie,bas_ie) = Matrices.S3_B_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.S4_B_EP_loc(bas_ie,bas_ie) = Matrices.S4_B_EP_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .*face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

        Matrices.BT1_EP_loc(bas_ie,bas_ie) = Matrices.BT1_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* face.nx .* face.gradedgeqx)' *  face.phiedgeq;
        Matrices.BT2_EP_loc(bas_ie,bas_ie) = Matrices.BT2_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* face.ny .* face.gradedgeqx)' *  face.phiedgeq;
        Matrices.BT3_EP_loc(bas_ie,bas_ie) = Matrices.BT3_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* face.nx .* face.gradedgeqy)' *  face.phiedgeq;
        Matrices.BT4_EP_loc(bas_ie,bas_ie) = Matrices.BT4_EP_loc(bas_ie,bas_ie) + (face.ds .* m_ave .* face.ny .* face.gradedgeqy)' *  face.phiedgeq;

    %% Internal faces with acoustic neighbor
    elseif face.neigh_ie > 0 && femregion.label(face.neigh_ie) == 'A'
    
        %% Element itself
        Matrices.D1_loc(bas_ie,bas_ie) = Matrices.D1_loc(bas_ie,bas_ie) + (face.ds .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.D2_loc(bas_ie,bas_ie) = Matrices.D2_loc(bas_ie,bas_ie) + (face.ds .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
        Matrices.D3_loc(bas_ie,bas_ie) = Matrices.D3_loc(bas_ie,bas_ie) + (face.ds .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
        Matrices.D4_loc(bas_ie,bas_ie) = Matrices.D4_loc(bas_ie,bas_ie) + (face.ds .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;


        if (Data.tau == 0)

            Matrices.S1_B_loc(bas_ie,bas_ie) = Matrices.S1_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S2_B_loc(bas_ie,bas_ie) = Matrices.S2_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S3_B_loc(bas_ie,bas_ie) = Matrices.S3_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S4_B_loc(bas_ie,bas_ie) = Matrices.S4_B_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

            Matrices.S1_B_beta_loc(bas_ie,bas_ie) = Matrices.S1_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S2_B_beta_loc(bas_ie,bas_ie) = Matrices.S2_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S3_B_beta_loc(bas_ie,bas_ie) = Matrices.S3_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S4_B_beta_loc(bas_ie,bas_ie) = Matrices.S4_B_beta_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

            Matrices.S1_B_beta2_loc(bas_ie,bas_ie) = Matrices.S1_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.nx .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S2_B_beta2_loc(bas_ie,bas_ie) = Matrices.S2_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.ny .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S3_B_beta2_loc(bas_ie,bas_ie) = Matrices.S3_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.nx .* face.phiedgeq)' * face.phiedgeq;
            Matrices.S4_B_beta2_loc(bas_ie,bas_ie) = Matrices.S4_B_beta2_loc(bas_ie,bas_ie) + face.penalty_geom.max(iedg) * (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.ny .* face.phiedgeq)' * face.phiedgeq;

            Matrices.BT1_loc(bas_ie,bas_ie) = Matrices.BT1_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* face.nx .* face.gradedgeqx)' * face.phiedgeq;
            Matrices.BT2_loc(bas_ie,bas_ie) = Matrices.BT2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* face.ny .* face.gradedgeqx)' * face.phiedgeq;
            Matrices.BT3_loc(bas_ie,bas_ie) = Matrices.BT3_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* face.nx .* face.gradedgeqy)' * face.phiedgeq;
            Matrices.BT4_loc(bas_ie,bas_ie) = Matrices.BT4_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* face.ny .* face.gradedgeqy)' * face.phiedgeq;

            Matrices.BT1_beta_loc(bas_ie,bas_ie) = Matrices.BT1_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave .* face.nx .* face.gradedgeqx)' * face.phiedgeq;
            Matrices.BT2_beta_loc(bas_ie,bas_ie) = Matrices.BT2_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave .* face.ny .* face.gradedgeqx)' * face.phiedgeq;
            Matrices.BT3_beta_loc(bas_ie,bas_ie) = Matrices.BT3_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave .* face.nx .* face.gradedgeqy)' * face.phiedgeq;
            Matrices.BT4_beta_loc(bas_ie,bas_ie) = Matrices.BT4_beta_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave .* face.ny .* face.gradedgeqy)' * face.phiedgeq;

            Matrices.BT1_beta2_loc(bas_ie,bas_ie) = Matrices.BT1_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.gradedgeqx)' * face.phiedgeq;
            Matrices.BT2_beta2_loc(bas_ie,bas_ie) = Matrices.BT2_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.gradedgeqx)' * face.phiedgeq;
            Matrices.BT3_beta2_loc(bas_ie,bas_ie) = Matrices.BT3_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave.^2 .* face.nx .* face.gradedgeqy)' * face.phiedgeq;
            Matrices.BT4_beta2_loc(bas_ie,bas_ie) = Matrices.BT4_beta2_loc(bas_ie,bas_ie) + 0.5 * (face.ds .* m_ave .* beta_ave.^2 .* face.ny .* face.gradedgeqy)' * face.phiedgeq;

        end

    end
     
end
