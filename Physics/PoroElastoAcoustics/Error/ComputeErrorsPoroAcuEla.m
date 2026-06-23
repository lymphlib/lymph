%> @file ComputeErrorsPoroAcuEla.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Compute errors for convergence analysis
%>
%==========================================================================
%> @section classComputeErrorsPoroAcuEla description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data        Struct with problem's data
%> @param mesh	      Struct containing mesh information (region+neighbor)	
%> @param femregion   Struct containing all the information 
%>                    about the finite element approximation
%> @param Solution    Problem's solution structure
%> @param time        Current time
%>
%> @retval Error      Structure with computed errors 
%>                   
%==========================================================================

function [Error] = ComputeErrorsPoroAcuEla(Data, mesh, femregion, Solution, time)
    
    Funcs.Preallocation    = @ErrorPreallocationElastodynamics;
    Funcs.VolumeAssemblyST = @VolumeErrorAssemblyElastodynamics;
    Funcs.FacesAssembly    = @FacesErrorAssemblyElastodynamics; 
    Funcs.FinalMatrices    = @ErrorElastodynamics;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;

    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = false;

    AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
    AssembInfo.ass_face_vec = ones(femregion.nel,1);
    
    AssembInfo.Solutions = Solution;
    AssembInfo.t = time;
    AssembInfo.label = femregion.label;
    AssembInfo.nbases = femregion.nbases;

    [Error] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Error] = ErrorPreallocationElastodynamics(GenMatrices)
    Error.Volume.err_ue_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_up_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_wp_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_phi_L2_loc  = GenMatrices.CellVector;
    Error.Volume.err_p_L2_loc    = GenMatrices.CellVector;

    Error.Volume.err_B_loc       = GenMatrices.CellVector;

    Error.Volume.err_dot_ue_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_dot_up_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_dot_wp_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_dot_uwp_L2_loc  = GenMatrices.CellVector;
    Error.Volume.err_dot_phi_L2_loc  = GenMatrices.CellVector;
    
    Error.Volume.err_dGep_loc     = GenMatrices.CellVector;
    Error.Volume.err_dGp_w_loc    = GenMatrices.CellVector;
    Error.Volume.err_dGp_beta_loc = GenMatrices.CellVector;
    Error.Volume.err_dGe_loc      = GenMatrices.CellVector;
    Error.Volume.err_dGa_loc      = GenMatrices.CellVector;

    Error.Faces.err_dGep_jumps_loc      = GenMatrices.CellVector;
    Error.Faces.err_dGp_w_jumps_loc     = GenMatrices.CellVector;
    Error.Faces.err_dGp_beta_jumps_loc  = GenMatrices.CellVector;
    Error.Faces.err_dGe_jumps_loc       = GenMatrices.CellVector;
    Error.Faces.err_dGa_jumps_loc       = GenMatrices.CellVector;

    Error.Faces.err_interf_PA_loc       = GenMatrices.CellVector;
end


%% Error assembly functions
function [Error] = VolumeErrorAssemblyElastodynamics(Data, Error, elem, ie, id, ~, AssembInfo)
         
    idx_in = AssembInfo.nbases(1:ie-1)'*(AssembInfo.label(1:ie-1) == AssembInfo.label(ie))+1;
    idx_end = idx_in+AssembInfo.nbases(ie)-1;
    elem.index = idx_in:idx_end;
    if AssembInfo.label(ie) == 'P' 
        
        % Evaluation of physical parameters
        lam     =  Data.lam{id}(elem.xq,elem.yq);
        mu       = Data.mu{id}(elem.xq,elem.yq);
        rho      = Data.rho{id}(elem.xq,elem.yq);
        rho_w    = Data.rho_w{id}(elem.xq,elem.yq);
        rho_f    = Data.rho_f{id}(elem.xq,elem.yq);
        m        = Data.m{id}(elem.xq,elem.yq);
        beta     = Data.beta{id}(elem.xq,elem.yq);
        eta_kper = Data.eta{id}(elem.xq,elem.yq)./Data.k_per{id}(elem.xq,elem.yq);
        
        % Exact solution evaluation
        up1ex_loc = Data.up_ex{1}(elem.xq,elem.yq).*Data.up_t_ex{1}(AssembInfo.t);
        up2ex_loc = Data.up_ex{2}(elem.xq,elem.yq).*Data.up_t_ex{1}(AssembInfo.t);
        dot_up1ex_loc = Data.up_ex{2}(elem.xq,elem.yq).*Data.dup_t_ex{1}(AssembInfo.t);
        dot_up2ex_loc = Data.up_ex{2}(elem.xq,elem.yq).*Data.dup_t_ex{1}(AssembInfo.t);
 
        gradx_up1ex_loc = Data.grad_up_ex{1}(elem.xq,elem.yq).*Data.up_t_ex{1}(AssembInfo.t);
        grady_up1ex_loc = Data.grad_up_ex{2}(elem.xq,elem.yq).*Data.up_t_ex{1}(AssembInfo.t);
        gradx_up2ex_loc = Data.grad_up_ex{3}(elem.xq,elem.yq).*Data.up_t_ex{1}(AssembInfo.t);
        grady_up2ex_loc = Data.grad_up_ex{4}(elem.xq,elem.yq).*Data.up_t_ex{1}(AssembInfo.t);

        wp1ex_loc = Data.wp_ex{1}(elem.xq,elem.yq).*Data.wp_t_ex{1}(AssembInfo.t);
        wp2ex_loc = Data.wp_ex{2}(elem.xq,elem.yq).*Data.wp_t_ex{1}(AssembInfo.t);
        dot_wp1ex_loc = Data.wp_ex{2}(elem.xq,elem.yq).*Data.dwp_t_ex{1}(AssembInfo.t);
        dot_wp2ex_loc = Data.wp_ex{2}(elem.xq,elem.yq).*Data.dwp_t_ex{1}(AssembInfo.t);
 
        gradx_wp1ex_loc = Data.grad_wp_ex{1}(elem.xq,elem.yq).*Data.wp_t_ex{1}(AssembInfo.t);
        grady_wp1ex_loc = Data.grad_wp_ex{2}(elem.xq,elem.yq).*Data.wp_t_ex{1}(AssembInfo.t);
        gradx_wp2ex_loc = Data.grad_wp_ex{3}(elem.xq,elem.yq).*Data.wp_t_ex{1}(AssembInfo.t);
        grady_wp2ex_loc = Data.grad_wp_ex{4}(elem.xq,elem.yq).*Data.wp_t_ex{1}(AssembInfo.t);

        % Numerical solution evaluation
        up1h_loc     = elem.phiq * AssembInfo.Solutions.up_h(elem.index);
        up2h_loc     = elem.phiq * AssembInfo.Solutions.up_h(end/2+elem.index);
        dot_up1h_loc = elem.phiq * AssembInfo.Solutions.dot_up_h(elem.index);
        dot_up2h_loc = elem.phiq * AssembInfo.Solutions.dot_up_h(end/2+elem.index);
        
        gradx_up1h_loc = elem.gradqx * AssembInfo.Solutions.up_h(elem.index);
        grady_up1h_loc = elem.gradqy * AssembInfo.Solutions.up_h(end/2+elem.index);
        gradx_up2h_loc = elem.gradqx * AssembInfo.Solutions.up_h(elem.index);
        grady_up2h_loc = elem.gradqy * AssembInfo.Solutions.up_h(end/2+elem.index);
        
        wp1h_loc     = elem.phiq * AssembInfo.Solutions.wp_h(elem.index);
        wp2h_loc     = elem.phiq * AssembInfo.Solutions.wp_h(end/2+elem.index);
        dot_wp1h_loc = elem.phiq * AssembInfo.Solutions.dot_wp_h(elem.index);
        dot_wp2h_loc = elem.phiq * AssembInfo.Solutions.dot_wp_h(end/2+elem.index);
        
        gradx_wp1h_loc = elem.gradqx * AssembInfo.Solutions.wp_h(elem.index);
        grady_wp1h_loc = elem.gradqy * AssembInfo.Solutions.wp_h(end/2+elem.index);
        gradx_wp2h_loc = elem.gradqx * AssembInfo.Solutions.wp_h(elem.index);
        grady_wp2h_loc = elem.gradqy * AssembInfo.Solutions.wp_h(end/2+elem.index);
        
        ph_loc = m .* (gradx_wp1h_loc + grady_wp2h_loc + beta .* (gradx_up1h_loc + grady_up2h_loc));

        % Local error integral assembly
        Error.err_up_L2_loc = (elem.dx .* (up1h_loc - up1ex_loc))' * (up1h_loc - up1ex_loc) ...
                            + (elem.dx .* (up2h_loc - up2ex_loc))' * (up2h_loc - up2ex_loc); 
       
        Error.err_dot_up_L2_loc = (elem.dx .* rho .* (dot_up1h_loc - dot_up1ex_loc))' * (dot_up1h_loc - dot_up1ex_loc) ...
                                + (elem.dx .* rho .* (dot_up2h_loc - dot_up2ex_loc))' * (dot_up2h_loc - dot_up2ex_loc); 
        
        Error.err_B_loc = (elem.dx .* eta_kper .* (wp1h_loc - wp1ex_loc))' * (wp1h_loc - wp1ex_loc) ...
                        + (elem.dx .* eta_kper .* (wp2h_loc - wp2ex_loc))' * (wp2h_loc - wp2ex_loc); 
       
        Error.err_wp_L2_loc = (elem.dx .* (wp1h_loc - wp1ex_loc))' * (wp1h_loc - wp1ex_loc) ...
                            + (elem.dx .* (wp2h_loc - wp2ex_loc))' * (wp2h_loc - wp2ex_loc); 

        Error.err_dot_wp_L2_loc = (elem.dx .* rho_w .* (dot_wp1h_loc - dot_wp1ex_loc))' * (dot_wp1h_loc - dot_wp1ex_loc) ...
                                + (elem.dx .* rho_w .* (dot_wp2h_loc - dot_wp2ex_loc))' * (dot_wp2h_loc - dot_wp2ex_loc); 
        
        Error.err_dot_uwp_L2_loc = (elem.dx .* rho_f .* (dot_wp1h_loc - dot_wp1ex_loc))' * (dot_wp1h_loc - dot_wp1ex_loc) ...
                                 + (elem.dx .* rho_f .* (dot_wp2h_loc - dot_wp2ex_loc))' * (dot_wp2h_loc - dot_wp2ex_loc); 

        Error.err_p_L2_loc = (elem.dx .* ph_loc)' * ph_loc;

        Error.err_dGep_loc  = (elem.dx .* (lam+2*mu) .* (gradx_up1h_loc - gradx_up1ex_loc))' * (gradx_up1h_loc - gradx_up1ex_loc) ...
                            + (elem.dx .* mu  .* (grady_up1h_loc - grady_up1ex_loc))' * (grady_up1h_loc - grady_up1ex_loc) ...
                            + (elem.dx .* lam .* (gradx_up1h_loc - gradx_up1ex_loc))' * (grady_up2h_loc - grady_up2ex_loc) ...
                            + (elem.dx .* mu  .* (grady_up1h_loc - grady_up1ex_loc))' * (gradx_up2h_loc - gradx_up2ex_loc) ...
                            + (elem.dx .* lam .* (grady_up2h_loc - grady_up2ex_loc))' * (gradx_up1h_loc - gradx_up1ex_loc) ...
                            + (elem.dx .* mu  .* (gradx_up2h_loc - gradx_up2ex_loc))' * (grady_up1h_loc - grady_up1ex_loc) ...
                            + (elem.dx .* (lam+2*mu) .* (grady_up2h_loc - grady_up2ex_loc))' * (grady_up2h_loc - grady_up2ex_loc) ...
                            + (elem.dx .* mu  .* (gradx_up2h_loc - gradx_up2ex_loc))' * (gradx_up2h_loc - gradx_up2ex_loc);

        Error.err_dGp_w_loc = (elem.dx .* m .* (gradx_wp1h_loc - gradx_wp1ex_loc))' * (gradx_wp1h_loc - gradx_wp1ex_loc) ...
                            + (elem.dx .* m .* (grady_wp1h_loc - grady_wp1ex_loc))' * (grady_wp1h_loc - grady_wp1ex_loc) ...
                            + (elem.dx .* m .* (gradx_wp1h_loc - gradx_wp1ex_loc))' * (grady_wp2h_loc - grady_wp2ex_loc) ...
                            + (elem.dx .* m .* (grady_wp1h_loc - grady_wp1ex_loc))' * (gradx_wp2h_loc - gradx_wp2ex_loc) ...
                            + (elem.dx .* m .* (grady_wp2h_loc - grady_wp2ex_loc))' * (gradx_wp1h_loc - gradx_wp1ex_loc) ...
                            + (elem.dx .* m .* (gradx_wp2h_loc - gradx_wp2ex_loc))' * (grady_wp1h_loc - grady_wp1ex_loc) ...
                            + (elem.dx .* m .* (grady_wp2h_loc - grady_wp2ex_loc))' * (grady_wp2h_loc - grady_wp2ex_loc) ...
                            + (elem.dx .* m .* (gradx_wp2h_loc - gradx_wp2ex_loc))' * (gradx_wp2h_loc - gradx_wp2ex_loc);

        Error.err_dGp_beta_loc  = (elem.dx .* (m .* beta) .* (gradx_up1h_loc - gradx_up1ex_loc))' * (gradx_up1h_loc - gradx_up1ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (grady_up1h_loc - grady_up1ex_loc))' * (grady_up1h_loc - grady_up1ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (gradx_up1h_loc - gradx_up1ex_loc))' * (grady_up2h_loc - grady_up2ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (grady_up1h_loc - grady_up1ex_loc))' * (gradx_up2h_loc - gradx_up2ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (grady_up2h_loc - grady_up2ex_loc))' * (gradx_up1h_loc - gradx_up1ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (gradx_up2h_loc - gradx_up2ex_loc))' * (grady_up1h_loc - grady_up1ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (grady_up2h_loc - grady_up2ex_loc))' * (grady_up2h_loc - grady_up2ex_loc) ...
                                + (elem.dx .* (m .* beta) .* (gradx_up2h_loc - gradx_up2ex_loc))' * (gradx_up2h_loc - gradx_up2ex_loc);

    elseif AssembInfo.label(ie) == 'E'
        % Evaluation of physical parameters
        rho_el = Data.rho_el{id}(elem.xq,elem.yq);
        lam    = Data.lam_el{id}(elem.xq,elem.yq);
        mu     = Data.mu_el{id}(elem.xq,elem.yq);

        % Exact solution evaluation
        ue1ex_loc = Data.ue_ex{1}(elem.xq,elem.yq).*Data.ue_t_ex{1}(AssembInfo.t);
        ue2ex_loc = Data.ue_ex{2}(elem.xq,elem.yq).*Data.ue_t_ex{1}(AssembInfo.t);
        dot_ue1ex_loc = Data.ue_ex{2}(elem.xq,elem.yq).*Data.due_t_ex{1}(AssembInfo.t);
        dot_ue2ex_loc = Data.ue_ex{2}(elem.xq,elem.yq).*Data.due_t_ex{1}(AssembInfo.t);
        
        gradx_ue1ex_loc = Data.grad_ue_ex{1}(elem.xq,elem.yq).*Data.ue_t_ex{1}(AssembInfo.t);
        grady_ue1ex_loc = Data.grad_ue_ex{2}(elem.xq,elem.yq).*Data.ue_t_ex{1}(AssembInfo.t);
        gradx_ue2ex_loc = Data.grad_ue_ex{3}(elem.xq,elem.yq).*Data.ue_t_ex{1}(AssembInfo.t);
        grady_ue2ex_loc = Data.grad_ue_ex{4}(elem.xq,elem.yq).*Data.ue_t_ex{1}(AssembInfo.t);
    
        % Numerical solution evaluation
        ue1h_loc     = elem.phiq * AssembInfo.Solutions.ue_h(elem.index);
        ue2h_loc     = elem.phiq * AssembInfo.Solutions.ue_h(end/2+elem.index);
        dot_ue1h_loc = elem.phiq * AssembInfo.Solutions.dot_ue_h(elem.index);
        dot_ue2h_loc = elem.phiq * AssembInfo.Solutions.dot_ue_h(end/2+elem.index);
        
        gradx_ue1h_loc = elem.gradqx * AssembInfo.Solutions.ue_h(elem.index);
        grady_ue1h_loc = elem.gradqy * AssembInfo.Solutions.ue_h(end/2+elem.index);
        gradx_ue2h_loc = elem.gradqx * AssembInfo.Solutions.ue_h(elem.index);
        grady_ue2h_loc = elem.gradqy * AssembInfo.Solutions.ue_h(end/2+elem.index);
        
        % Local error integral assembly
        Error.err_ue_L2_loc = (elem.dx .* (ue1h_loc - ue1ex_loc))' * (ue1h_loc - ue1ex_loc) ...
                            + (elem.dx .* (ue2h_loc - ue2ex_loc))' * (ue2h_loc - ue2ex_loc); 
       
        Error.err_dot_ue_L2_loc = (elem.dx .* rho_el .* (dot_ue1h_loc - dot_ue1ex_loc))' * (dot_ue1h_loc - dot_ue1ex_loc) ...
                                + (elem.dx .* rho_el .* (dot_ue2h_loc - dot_ue2ex_loc))' * (dot_ue2h_loc - dot_ue2ex_loc); 
              
        Error.err_dGe_loc  = (elem.dx .* (lam+2*mu) .* (gradx_ue1h_loc - gradx_ue1ex_loc))' * (gradx_ue1h_loc - gradx_ue1ex_loc) ...
                           + (elem.dx .* mu  .* (grady_ue1h_loc - grady_ue1ex_loc))' * (grady_ue1h_loc - grady_ue1ex_loc) ...
                           + (elem.dx .* lam .* (gradx_ue1h_loc - gradx_ue1ex_loc))' * (grady_ue2h_loc - grady_ue2ex_loc) ...
                           + (elem.dx .* mu  .* (grady_ue1h_loc - grady_ue1ex_loc))' * (gradx_ue2h_loc - gradx_ue2ex_loc) ...
                           + (elem.dx .* lam .* (grady_ue2h_loc - grady_ue2ex_loc))' * (gradx_ue1h_loc - gradx_ue1ex_loc) ...
                           + (elem.dx .* mu  .* (gradx_ue2h_loc - gradx_ue2ex_loc))' * (grady_ue1h_loc - grady_ue1ex_loc) ...
                           + (elem.dx .* (lam+2*mu) .* (grady_ue2h_loc - grady_ue2ex_loc))' * (grady_ue2h_loc - grady_ue2ex_loc) ...
                           + (elem.dx .* mu  .* (gradx_ue2h_loc - gradx_ue2ex_loc))' * (gradx_ue2h_loc - gradx_ue2ex_loc);

    elseif AssembInfo.label(ie) == 'A'

        rho_a   = Data.rho_a{id}(elem.xq,elem.yq);
        c       = Data.c{id}(elem.xq,elem.yq);

        % Exact solution evaluation
        phiex_loc     = Data.phi_ex{1}(elem.xq,elem.yq).*Data.phi_t_ex{1}(AssembInfo.t);
        dot_phiex_loc = Data.phi_ex{1}(elem.xq,elem.yq).*Data.dphi_t_ex{1}(AssembInfo.t);
        gradx_phiex_loc  = Data.grad_phi_ex{1}(elem.xq,elem.yq).*Data.phi_t_ex{1}(AssembInfo.t);
        grady_phiex_loc  = Data.grad_phi_ex{2}(elem.xq,elem.yq).*Data.phi_t_ex{1}(AssembInfo.t);
        
        % Numerical solution evaluation
        phih_loc     = elem.phiq * AssembInfo.Solutions.phi_h(elem.index);
        dot_phih_loc = elem.phiq * AssembInfo.Solutions.dot_phi_h(elem.index);
        gradx_phih_loc  = elem.gradqx * AssembInfo.Solutions.phi_h(elem.index);
        grady_phih_loc  = elem.gradqy * AssembInfo.Solutions.phi_h(elem.index);
        
        % Local error integral assembly
        Error.err_phi_L2_loc = (elem.dx .* (phih_loc - phiex_loc))' * (phih_loc - phiex_loc);
        Error.err_dot_phi_L2_loc = (elem.dx .* (c.^(-2) .* rho_a) .* (dot_phih_loc - dot_phiex_loc))' * (dot_phih_loc - dot_phiex_loc);
        Error.err_dGa_loc = (elem.dx .* rho_a .* (gradx_phih_loc - gradx_phiex_loc))' * (gradx_phih_loc - gradx_phiex_loc) ...
                          + (elem.dx .* rho_a .* (grady_phih_loc - grady_phiex_loc))' * (grady_phih_loc - grady_phiex_loc);
        
    end
end

function [Error] = FacesErrorAssemblyElastodynamics(Data, femregion, Error, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;
    
    idx_in = AssembInfo.nbases(1:ie-1)'*(AssembInfo.label(1:ie-1) == AssembInfo.label(ie))+1;
    idx_end = idx_in+AssembInfo.nbases(ie)-1;
    index_self = idx_in:idx_end;
    
    if face.neigh_ie > 0
        if femregion.label(ie) == femregion.label(face.neigh_ie)
        idx_in = AssembInfo.nbases(1:face.neigh_ie-1)'*(AssembInfo.label(1:face.neigh_ie-1) == AssembInfo.label(face.neigh_ie))+1;
            idx_end = idx_in+AssembInfo.nbases(face.neigh_ie)-1;
            index_neigh = idx_in:idx_end;
        end
    end

    if AssembInfo.label(ie) == 'P'

        mu   = Data.mu{femregion.id_phys(ie)}(face.xq,face.yq);
        lam  = Data.lam{femregion.id_phys(ie)}(face.xq,face.yq);
        m    = Data.m{femregion.id_phys(ie)}(face.xq,face.yq);
        beta = Data.beta{femregion.id_phys(ie)}(face.xq,face.yq);
        
        %Evaluation solution
        up1h_self = face.phiedgeq * AssembInfo.Solutions.up_h(index_self);
        up2h_self = face.phiedgeq * AssembInfo.Solutions.up_h(end/2+index_self);
        wp1h_self = face.phiedgeq * AssembInfo.Solutions.wp_h(index_self);
        wp2h_self = face.phiedgeq * AssembInfo.Solutions.wp_h(end/2+index_self);
            
        %% Dirichlet boundary faces
        if face.neigh_ie == -1
    
            up1_ex = Data.up_ex{1}(face.xq,face.yq).*Data.up_t_ex{1}(AssembInfo.t);
            up2_ex = Data.up_ex{2}(face.xq,face.yq).*Data.up_t_ex{1}(AssembInfo.t);
            wp1_ex = Data.wp_ex{1}(face.xq,face.yq).*Data.wp_t_ex{1}(AssembInfo.t);
            wp2_ex = Data.wp_ex{2}(face.xq,face.yq).*Data.wp_t_ex{1}(AssembInfo.t);
            
            Error.err_dGep_jumps_loc = Error.err_dGep_jumps_loc ...
                                     + face.ds' * (face.penalty_geom.harm(iedg) * (lam + 2*mu) .*(up1h_self-up1_ex).^2) ...
                                     + face.ds' * (face.penalty_geom.harm(iedg) * (lam + 2*mu) .*(up2h_self-up2_ex).^2);
    
            Error.err_dGp_w_jumps_loc = Error.err_dGp_w_jumps_loc ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m .*(wp1h_self-wp1_ex).^2) ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m .*(wp2h_self-wp2_ex).^2);
           
            Error.err_dGp_beta_jumps_loc = Error.err_dGp_beta_jumps_loc ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m .* beta .*(up1h_self-up1_ex).^2) ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m .* beta .*(up2h_self-up2_ex).^2);
        
        %% Internal faces
        elseif face.neigh_ie > ie && AssembInfo.label(face.neigh_ie) == 'P'
    
            mu_n   = Data.mu{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
            lam_n  = Data.lam{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
            m_n    = Data.m{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
            beta_n = Data.beta{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
            
            lambda_ave = 2*lam .* lam_n ./ (lam + lam_n);
            mu_ave     = 2*mu .* mu_n ./ (mu + mu_n);
            m_ave      = 2*m .* m_n ./ (m + m_n);
            beta_ave   = 2*beta .* beta_n ./ (beta + beta_n);
             
            up1h_neigh = face.phiedgeqneigh * AssembInfo.Solutions.up_h(index_neigh);
            up2h_neigh = face.phiedgeqneigh * AssembInfo.Solutions.up_h(end/2+index_neigh);
            wp1h_neigh = face.phiedgeqneigh * AssembInfo.Solutions.wp_h(index_neigh);
            wp2h_neigh = face.phiedgeqneigh * AssembInfo.Solutions.wp_h(end/2+index_neigh);
             
            Error.err_dGep_jumps_loc = Error.err_dGe_jumps_loc ...
                                     + face.ds' * (face.penalty_geom.harm(iedg) * (lambda_ave + 2*mu_ave) .* (up1h_self-up1h_neigh).^2) ...
                                     + face.ds' * (face.penalty_geom.harm(iedg) * (lambda_ave + 2*mu_ave) .* (up2h_self-up2h_neigh).^2);
            
            Error.err_dGp_w_jumps_loc = Error.err_dGp_w_jumps_loc ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m_ave .* (wp1h_self-wp1h_neigh).^2) ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m_ave .* (wp2h_self-wp2h_neigh).^2);

            Error.err_dGp_beta_jumps_loc = Error.err_dGp_beta_jumps_loc ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m_ave .* beta_ave .* (up1h_self-up1h_neigh).^2) ...
                                      + face.ds' * (face.penalty_geom.harm(iedg) * m_ave .* beta_ave .* (up2h_self-up2h_neigh).^2);

        elseif face.neigh_ie > ie && AssembInfo.label(face.neigh_ie) == 'A'
    
            tau = Data.tau;
            
            if tau ~= 0
                wp1_ex = Data.wp_ex{1}(face.xq,face.yq).*Data.wp_t_ex{1}(AssembInfo.t);
                wp2_ex = Data.wp_ex{2}(face.xq,face.yq).*Data.wp_t_ex{1}(AssembInfo.t);
    
                Error.err_interf_PA_loc = Error.err_interf_PA_loc ...
                                         + face.ds' * ((1-tau)/tau .* ((wp1h_self-wp1_ex) * face.nx + (wp2h_self-wp2_ex) * face.ny).^2);
            end
        end

    elseif AssembInfo.label(ie) == 'E'
        
        mu  = Data.mu_el{femregion.id_phys(ie)}(face.xq,face.yq);
        lam = Data.lam_el{femregion.id_phys(ie)}(face.xq,face.yq);
    
        %Evaluation solution
        ue1h_self = face.phiedgeq * AssembInfo.Solutions.ue_h(index_self);
        ue2h_self = face.phiedgeq * AssembInfo.Solutions.ue_h(end/2+index_self);
            
        %% Dirichlet boundary faces
        if face.neigh_ie == -1
    
            ue1_ex = Data.ue_ex{1}(face.xq,face.yq).*Data.ue_t_ex{1}(AssembInfo.t);
            ue2_ex = Data.ue_ex{2}(face.xq,face.yq).*Data.ue_t_ex{1}(AssembInfo.t);
            
            Error.err_dGe_jumps_loc = Error.err_dGe_jumps_loc ...
                                   + face.ds' * (face.penalty_geom.harm(iedg) * (lam + 2*mu) .*(ue1h_self-ue1_ex).^2) ...
                                   + face.ds' * (face.penalty_geom.harm(iedg) * (lam + 2*mu) .*(ue2h_self-ue2_ex).^2);
    
        %% Internal faces
        elseif face.neigh_ie > ie && AssembInfo.label(face.neigh_ie) == 'E'
    
            mu_n  = Data.mu_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
            lam_n = Data.lam_el{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
     
            lambda_ave = 2*lam .* lam_n ./ (lam + lam_n);
            mu_ave     = 2*mu .* mu_n ./ (mu + mu_n);
             
            ue1h_neigh = face.phiedgeqneigh * AssembInfo.Solutions.ue_h(index_neigh);
            ue2h_neigh = face.phiedgeqneigh * AssembInfo.Solutions.ue_h(end/2+index_neigh);
             
            Error.err_dGe_jumps_loc = Error.err_dGe_jumps_loc ...
                                    + face.ds' * (face.penalty_geom.harm(iedg) * (lambda_ave + 2*mu_ave) .* (ue1h_self-ue1h_neigh).^2) ...
                                    + face.ds' * (face.penalty_geom.harm(iedg) * (lambda_ave + 2*mu_ave) .* (ue2h_self-ue2h_neigh).^2);
        
        end

    elseif AssembInfo.label(ie) == 'A'
        
        
        rho_a = Data.rho_a{femregion.id_phys(ie)}(face.xq,face.yq);
 
        %Evaluation solution
        phih_self = face.phiedgeq * AssembInfo.Solutions.phi_h(index_self);
       
        %% Dirichlet boundary faces
        if face.neigh_ie == -1

            phiex_loc = Data.phi_ex{1}(face.xq,face.yq).*Data.phi_t_ex{1}(AssembInfo.t);
            
            Error.err_dGa_jumps_loc = Error.err_dGa_jumps_loc ...
                                    + face.ds' * (face.penalty_geom.harm(iedg) * rho_a .*(phih_self-phiex_loc).^2);

            %% Internal faces
        elseif face.neigh_ie > ie && AssembInfo.label(face.neigh_ie) == 'A'

            rho_a_n   = Data.rho_a{femregion.id_phys(ie)}(face.xq,face.yq);
            rho_a_ave = 2*rho_a .* rho_a_n ./ (rho_a + rho_a_n);
            
            phih_neigh = face.phiedgeqneigh * AssembInfo.Solutions.phi_h(index_neigh);
            
            Error.err_dGa_jumps_loc = Error.err_dGa_jumps_loc ...
                                    + face.ds' * (face.penalty_geom.harm(iedg) * rho_a_ave .* (phih_self-phih_neigh).^2);
        end
    end

end

%% Final error construction function
function [Error] = ErrorElastodynamics(Error_loc)
    Error.error_L2_ue  = sqrt(sum(Error_loc.Volume.err_ue_L2_loc));
    Error.error_L2_up  = sqrt(sum(Error_loc.Volume.err_up_L2_loc));
    Error.error_L2_wp  = sqrt(sum(Error_loc.Volume.err_wp_L2_loc));
    Error.error_L2_phi = sqrt(sum(Error_loc.Volume.err_phi_L2_loc));
    
    Error.error_L2_d        = sqrt(Error.error_L2_ue^2 + Error.error_L2_up^2 + Error.error_L2_wp^2 + Error.error_L2_phi^2);
    Error.error_L2_p        = sqrt(sum(Error_loc.Volume.err_p_L2_loc));
    Error.error_L2_pressure = sqrt(sum(Error_loc.Volume.err_p_L2_loc + Error_loc.Volume.err_dot_phi_L2_loc));

    Error.error_L2_dot_ue  = sqrt(sum(Error_loc.Volume.err_dot_ue_L2_loc));
    Error.error_L2_dot_up  = sqrt(sum(Error_loc.Volume.err_dot_up_L2_loc));
    Error.error_L2_dot_wp  = sqrt(sum(Error_loc.Volume.err_dot_wp_L2_loc));
    Error.error_L2_dot_uwp = sqrt(sum(Error_loc.Volume.err_dot_uwp_L2_loc));
    Error.error_L2_dot_phi = sqrt(sum(Error_loc.Volume.err_dot_phi_L2_loc));
    
    Error.error_L2_v = sqrt(Error.error_L2_dot_ue^2 + Error.error_L2_dot_up^2 + Error.error_L2_dot_wp^2 + Error.error_L2_dot_uwp^2 + Error.error_L2_dot_phi^2);
    
    Error.error_dGep  = sqrt(sum(Error_loc.Volume.err_dGep_loc+Error_loc.Faces.err_dGe_jumps_loc));
    Error.error_dGp_w = sqrt(sum(Error_loc.Volume.err_dGp_w_loc+Error_loc.Faces.err_dGp_w_jumps_loc));
    Error.error_dGp   = sqrt(sum(Error_loc.Volume.err_dGp_w_loc+Error_loc.Faces.err_dGp_w_jumps_loc ... 
                               + Error_loc.Volume.err_dGp_beta_loc+Error_loc.Faces.err_dGp_beta_jumps_loc));
    Error.error_dGe   = sqrt(sum(Error_loc.Volume.err_dGe_loc+Error_loc.Faces.err_dGe_jumps_loc));
    Error.error_dGa   = sqrt(sum(Error_loc.Volume.err_dGa_loc+Error_loc.Faces.err_dGa_jumps_loc));
    
    Error.error_dG    = sqrt(Error.error_dGep^2 + Error.error_dGp^2 + Error.error_dGa^2 + Error.error_dGe^2);
    
    Error.error_B_interface = sqrt(sum(Error_loc.Faces.err_interf_PA_loc));
    Error.error_B = sqrt(sum(Error_loc.Volume.err_B_loc+Error_loc.Faces.err_interf_PA_loc));

    Error.error_Energy = sqrt(Error.error_L2_v^2 + Error.error_dG^2 + Error.error_B^2);

end
