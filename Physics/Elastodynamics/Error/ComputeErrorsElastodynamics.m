%> @file ComputeErrorsElastodynamics.m
%> @author Mattia Corti
%> @date 8 May 2026
%> @brief Compute errors for convergence analysis
%>
%==========================================================================
%> @section classComputeErrorsElastodynamics description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data        Struct with problem's data
%> @param mesh	      Struct containing mesh information (region+neighbor)	
%> @param femregion   Struct containing all the information 
%>                    about the finite element approximation
%> @param U           Problem's solution
%> @param time        Current time
%>
%> @retval Error      Structure with computed L2 and dG errors 
%>                   
%==========================================================================

function [Error] = ComputeErrorsElastodynamics(Data, mesh, femregion, U, time)
    
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
    
    AssembInfo.ue1_h = U(1:femregion.ndof);
    AssembInfo.ue2_h = U(femregion.ndof+1:2*femregion.ndof);
    AssembInfo.ve1_h = U(2*femregion.ndof+1:3*femregion.ndof);
    AssembInfo.ve2_h = U(3*femregion.ndof+1:end);
    AssembInfo.t = time;

    [Error] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Error] = ErrorPreallocationElastodynamics(GenMatrices)
    Error.Volume.err_u_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_v_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_u_H1_loc   = GenMatrices.CellVector;

    Error.Faces.err_u_jumps_loc = GenMatrices.CellVector;
end


%% Error assembly functions
function [Error] = VolumeErrorAssemblyElastodynamics(Data, Error, elem, ie, id, ~, AssembInfo)
                
        % Evaluation of physical parameters
        rho_el = Data.rho_el{id}(elem.xq,elem.yq);
        lam    = Data.lam_el{id}(elem.xq,elem.yq);
        mu     = Data.mu_el{id}(elem.xq,elem.yq);

        % Exact solution evaluation
        ue1ex_loc = Data.ue_ex{1}(elem.xq,elem.yq,AssembInfo.t);
        ue2ex_loc = Data.ue_ex{2}(elem.xq,elem.yq,AssembInfo.t);
        ve1ex_loc = Data.due_dt_ex{1}(elem.xq,elem.yq,AssembInfo.t);
        ve2ex_loc = Data.due_dt_ex{2}(elem.xq,elem.yq,AssembInfo.t);
        gradx_u1ex_loc = Data.grad_ue_ex{1}(elem.xq,elem.yq, AssembInfo.t);
        grady_u1ex_loc = Data.grad_ue_ex{2}(elem.xq,elem.yq, AssembInfo.t);
        gradx_u2ex_loc = Data.grad_ue_ex{3}(elem.xq,elem.yq, AssembInfo.t);
        grady_u2ex_loc = Data.grad_ue_ex{4}(elem.xq,elem.yq, AssembInfo.t);
    
        % Numerical solution evaluation
        ue1h_loc = elem.phiq * AssembInfo.ue1_h(elem.index);
        ue2h_loc = elem.phiq * AssembInfo.ue2_h(elem.index);
        ve1h_loc = elem.phiq * AssembInfo.ve1_h(elem.index);
        ve2h_loc = elem.phiq * AssembInfo.ve2_h(elem.index);
        
        gradx_ue1h_loc = elem.gradqx * AssembInfo.ue1_h(elem.index);
        grady_ue1h_loc = elem.gradqy * AssembInfo.ue1_h(elem.index);
        gradx_ue2h_loc = elem.gradqx * AssembInfo.ue2_h(elem.index);
        grady_ue2h_loc = elem.gradqy * AssembInfo.ue2_h(elem.index);
        
        % Local error integral assembly
        Error.err_u_L2_loc = (elem.dx .* (ue1h_loc - ue1ex_loc))' * (ue1h_loc - ue1ex_loc) ...
                           + (elem.dx .* (ue2h_loc - ue2ex_loc))' * (ue2h_loc - ue2ex_loc); 
       
        Error.err_v_L2_loc = (elem.dx .* rho_el .* (ve1h_loc - ve1ex_loc))' * (ve1h_loc - ve1ex_loc) ...
                           + (elem.dx .* rho_el .* (ve2h_loc - ve2ex_loc))' * (ve2h_loc - ve2ex_loc); 
              
        Error.err_H1_loc   = (elem.dx .* (lam+2*mu) .* (gradx_ue1h_loc - gradx_u1ex_loc))' * (gradx_ue1h_loc - gradx_u1ex_loc) ...
                           + (elem.dx .* mu  .* (grady_ue1h_loc - grady_u1ex_loc))' * (grady_ue1h_loc - grady_u1ex_loc) ...
                           + (elem.dx .* lam .* (gradx_ue1h_loc - gradx_u1ex_loc))' * (grady_ue2h_loc - grady_u2ex_loc) ...
                           + (elem.dx .* mu  .* (grady_ue1h_loc - grady_u1ex_loc))' * (gradx_ue2h_loc - gradx_u2ex_loc) ...
                           + (elem.dx .* lam .* (grady_ue2h_loc - grady_u2ex_loc))' * (gradx_ue1h_loc - gradx_u1ex_loc) ...
                           + (elem.dx .* mu  .* (gradx_ue2h_loc - gradx_u2ex_loc))' * (grady_ue1h_loc - grady_u1ex_loc) ...
                           + (elem.dx .* (lam+2*mu) .* (grady_ue2h_loc - grady_u2ex_loc))' * (grady_ue2h_loc - grady_u2ex_loc) ...
                           + (elem.dx .* mu  .* (gradx_ue2h_loc - gradx_u2ex_loc))' * (gradx_ue2h_loc - gradx_u2ex_loc);

end

function [Error] = FacesErrorAssemblyElastodynamics(Data, femregion, Error, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;
    
    index_self  = sum(femregion.nbases(1:face.ie-1))+1:sum(femregion.nbases(1:face.ie));
    index_neigh = sum(femregion.nbases(1:face.neigh_ie-1))+1:sum(femregion.nbases(1:face.neigh_ie));

    mu    = Data.mu_el{femregion.id(ie)}(face.xq,face.yq);
    lam   = Data.lam_el{femregion.id(ie)}(face.xq,face.yq);

    %Evaluation solution
    ue1h_self = face.phiedgeq * AssembInfo.ue1_h(index_self);
    ue2h_self = face.phiedgeq * AssembInfo.ue2_h(index_self);
        
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        gD1 = Data.DirBCEla{1}(face.xq,face.yq,AssembInfo.t);
        gD2 = Data.DirBCEla{2}(face.xq,face.yq,AssembInfo.t);
        
        Error.err_u_jumps_loc = Error.err_u_jumps_loc ...
                              + face.ds' * (face.penalty_geom.harm(iedg) * (lam + 2*mu) .*(ue1h_self-gD1).^2) ...
                              + face.ds' * (face.penalty_geom.harm(iedg) * (lam + 2*mu) .*(ue2h_self-gD2).^2);

    %% Internal faces
    elseif face.neigh_ie > ie

        mu_n   = Data.mu_el{femregion.id(ie)}(face.xq,face.yq);
        lam_n  = Data.lam_el{femregion.id(ie)}(face.xq,face.yq);

        lambda_ave = 2*lam .* lam_n ./ (lam + lam_n);
        mu_ave     = 2*mu .* mu_n ./ (mu + mu_n);

        ue1h_neigh = face.phiedgeqneigh * AssembInfo.ue1_h(index_neigh);
        ue2h_neigh = face.phiedgeqneigh * AssembInfo.ue2_h(index_neigh);

        Error.err_u_jumps_loc = Error.err_u_jumps_loc ...
                              + face.ds' * (face.penalty_geom.harm(iedg) * (lambda_ave + 2*mu_ave) .* (ue1h_self-ue1h_neigh).^2) ...
                              + face.ds' * (face.penalty_geom.harm(iedg) * (lambda_ave + 2*mu_ave) .* (ue2h_self-ue2h_neigh).^2);

    end

end

%% Final error construction function
function [Error] = ErrorElastodynamics(Error_loc)
    Error.err_u_L2   = sqrt(sum(Error_loc.Volume.err_u_L2_loc));
    Error.err_u_dG   = sqrt(sum(Error_loc.Volume.err_H1_loc+Error_loc.Faces.err_u_jumps_loc));
    Error.err_v_L2   = sqrt(sum(Error_loc.Volume.err_v_L2_loc));
    Error.err_energy = sqrt(sum(Error_loc.Volume.err_v_L2_loc+Error_loc.Volume.err_H1_loc+Error_loc.Faces.err_u_jumps_loc));
end
