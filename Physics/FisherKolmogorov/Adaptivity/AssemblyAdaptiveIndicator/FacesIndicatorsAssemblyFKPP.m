%> @file   FacesIndicatorsAssemblyFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   18 February 2026
%> @brief Computation of the face indicators for FKPP equation.
%>
%==========================================================================
%> @section classFacesIndicatorsAssemblyFKPP Class description
%==========================================================================
%> @brief            Computation of the face indicators for FKPP equation.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Indicator  Indicator struct (to be modified) 
%>                    - \f$\tau_T\f$ = \f$\|h^1/2 \theta^1/2 [| \nabla u_h . t|] \|^2_{L2} + \| \sigma^1/2 \theta^1/2 \nabla (u_h - g) \|^2_{L2}  + 
%>                                     \|h^1/2 (1-\theta)^1/2 [| \nabla u_h^{old} . t|] \|^2_{L2} + \| \sigma^1/2 (1-theta)^1/2 \nabla (u_h^{old} - g) \|^2_{L2} \f$
%>                    - \f$\tau_J\f$ = \f$\|\sigma^1/2 \theta^1/2 [|u_h|] ||^2_{L2} + || \sigma^1/2 \theta^1/2 (u_h - g) ||^2_{L2} + 
%>                                     \|\sigma^1/2 (1-\theta)^1/2 [|u_h^{old}|] ||^2_{L2} + || sigma^1/2 (1-theta)^1/2 (u_h^{old} - g) ||^2_{L2}\f$
%>                    - \f$\tau_N\f$ = \f$\|h^1/2 \theta^1/2 [| \nabla u_h . n|] ||^2_{L2} + ||\sigma^1/2 \theta^1/2 \nabla (u_h - g_N)||^2_{L2} + 
%>                                     \|h^1/2 (1-\theta)^1/2 [| \nabla u_h^{old} . t|] ||^2_{L2} + || \sigma^1/2 (1-\theta)^1/2 \nabla (u_h^{old} - g_N) ||^2_{L2}\f$
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Indicator  Indicator struct
%                   
%==========================================================================

function [Indicator] = FacesIndicatorsAssemblyFKPP(Data, femregion, Indicator, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;
    
    index_self  = sum(femregion.nbases(1:face.ie-1))+1:sum(femregion.nbases(1:face.ie));
    index_neigh = sum(femregion.nbases(1:face.neigh_ie-1))+1:sum(femregion.nbases(1:face.neigh_ie));

    % Evaluation solution and gradients
    c_h_self  = face.phiedgeq * AssembInfo.c_h(index_self);
    grad_c_h_n_self  =  (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy) * AssembInfo.c_h(index_self);
    grad_c_h_t_self  = (-face.ny * face.gradedgeqx + face.nx * face.gradedgeqy) * AssembInfo.c_h(index_self);

    c_old_self  = face.phiedgeq * AssembInfo.c_old(index_self);
    grad_c_old_n_self  =  (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy) * AssembInfo.c_old(index_self);
    grad_c_old_t_self  = (-face.ny * face.gradedgeqx + face.nx * face.gradedgeqy) * AssembInfo.c_old(index_self);

    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        D_ext  = Data.D_ext{1}(face.xq,face.yq);
        
        gD  = Data.DirBC{1}(face.xq,face.yq,AssembInfo.t);
        gDx = Data.gradDirBC{1}(face.xq,face.yq,AssembInfo.t);
        gDy = Data.gradDirBC{2}(face.xq,face.yq,AssembInfo.t);

        gDold  = Data.DirBC{1}(face.xq,face.yq,AssembInfo.t-Data.dt);
        gDxold = Data.gradDirBC{1}(face.xq,face.yq,AssembInfo.t-Data.dt);
        gDyold = Data.gradDirBC{2}(face.xq,face.yq,AssembInfo.t-Data.dt);

        Indicator.tau_J = Indicator.tau_J + face.ds' * (face.penalty_geom.harm(iedg) * D_ext .* (Data.theta*(c_h_self - gD) + (1-Data.theta)*(c_old_self - gDold)).^2);
        Indicator.tau_T = Indicator.tau_T + face.ds' * (Data.h * (Data.theta*(grad_c_h_t_self - (-face.ny*gDx + face.nx*gDy)) + (1-Data.theta)*(grad_c_old_t_self - (-face.ny*gDxold + face.nx*gDyold))).^2);
   
    %% Internal faces
    elseif face.neigh_ie > 0

        D_ext  = Data.D_ext{1}(face.xq,face.yq);

        c_h_neigh        = face.phiedgeqneigh * AssembInfo.c_h(index_neigh);
        grad_c_h_n_neigh =  (face.nx * face.gradedgeqxneigh + face.ny * face.gradedgeqyneigh) * AssembInfo.c_h(index_neigh);
        grad_c_h_t_neigh = (-face.ny * face.gradedgeqxneigh + face.nx * face.gradedgeqyneigh) * AssembInfo.c_h(index_neigh);

        c_old_neigh        = face.phiedgeqneigh * AssembInfo.c_old(index_neigh);
        grad_c_old_n_neigh =  (face.nx * face.gradedgeqxneigh + face.ny * face.gradedgeqyneigh) * AssembInfo.c_old(index_neigh);
        grad_c_old_t_neigh = (-face.ny * face.gradedgeqxneigh + face.nx * face.gradedgeqyneigh) * AssembInfo.c_old(index_neigh);

        
        Indicator.tau_J = Indicator.tau_J + face.ds' * (face.penalty_geom.harm(iedg) * D_ext.* (Data.theta*(c_h_self - c_h_neigh) + (1-Data.theta)*(c_old_self - c_old_neigh)).^2);
        Indicator.tau_N = Indicator.tau_N + face.ds' * (Data.h * D_ext .* (Data.theta*(grad_c_h_n_self - grad_c_h_n_neigh) + (1-Data.theta)*(grad_c_old_n_self - grad_c_old_n_neigh)).^2);
        Indicator.tau_T = Indicator.tau_T + face.ds' * (Data.h * (Data.theta*(grad_c_h_t_self - grad_c_h_t_neigh) + (1-Data.theta)*(grad_c_old_t_self - grad_c_old_t_neigh)).^2);

    end

end
