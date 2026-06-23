%> @file   FacesIndicatorsAssemblyHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   18 February 2026
%> @brief Computation of the face indicators for heat equation.
%>
%==========================================================================
%> @section classFacesIndicatorsAssemblyHeat Class description
%==========================================================================
%> @brief            Computation of the face indicators for heat equation.
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

function [Indicator] = FacesIndicatorsAssemblyHeat(Data, femregion, Indicator, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;
    
    index_self  = sum(femregion.nbases(1:face.ie-1))+1:sum(femregion.nbases(1:face.ie));
    index_neigh = sum(femregion.nbases(1:face.neigh_ie-1))+1:sum(femregion.nbases(1:face.neigh_ie));

    % Evaluation solution and gradients
    uh_self  = face.phiedgeq * AssembInfo.uh(index_self);
    grad_uh_n_self  =  (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy) * AssembInfo.uh(index_self);
    grad_uh_t_self  = (-face.ny * face.gradedgeqx + face.nx * face.gradedgeqy) * AssembInfo.uh(index_self);

    uold_self  = face.phiedgeq * AssembInfo.u_old(index_self);
    grad_uold_n_self  =  (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy) * AssembInfo.u_old(index_self);
    grad_uold_t_self  = (-face.ny * face.gradedgeqx + face.nx * face.gradedgeqy) * AssembInfo.u_old(index_self);

    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        gD  = Data.DirBC{1}(face.xq,face.yq,AssembInfo.t);
        gDx = Data.gradDirBC{1}(face.xq,face.yq,AssembInfo.t);
        gDy = Data.gradDirBC{2}(face.xq,face.yq,AssembInfo.t);

        gDold  = Data.DirBC{1}(face.xq,face.yq,AssembInfo.t-Data.dt);
        gDxold = Data.gradDirBC{1}(face.xq,face.yq,AssembInfo.t-Data.dt);
        gDyold = Data.gradDirBC{2}(face.xq,face.yq,AssembInfo.t-Data.dt);

        Indicator.tau_J = Indicator.tau_J + face.ds' * (face.penalty_geom.max(iedg) * (Data.theta*(uh_self - gD) + (1-Data.theta)*(uold_self - gDold)).^2);
        Indicator.tau_T = Indicator.tau_T + face.ds' * (Data.h * (Data.theta*(grad_uh_t_self - (-face.ny*gDx + face.nx*gDy)) + (1-Data.theta)*(grad_uold_t_self - (-face.ny*gDxold + face.nx*gDyold))).^2);
   
    %% Internal faces
    elseif face.neigh_ie > 0

        uh_neigh        = face.phiedgeqneigh * AssembInfo.uh(index_neigh);
        grad_uh_n_neigh =  (face.nx * face.gradedgeqxneigh + face.ny * face.gradedgeqyneigh) * AssembInfo.uh(index_neigh);
        grad_uh_t_neigh = (-face.ny * face.gradedgeqxneigh + face.nx * face.gradedgeqyneigh) * AssembInfo.uh(index_neigh);

        uold_neigh        = face.phiedgeqneigh * AssembInfo.u_old(index_neigh);
        grad_uold_n_neigh =  (face.nx * face.gradedgeqxneigh + face.ny * face.gradedgeqyneigh) * AssembInfo.u_old(index_neigh);
        grad_uold_t_neigh = (-face.ny * face.gradedgeqxneigh + face.nx * face.gradedgeqyneigh) * AssembInfo.u_old(index_neigh);

        
        Indicator.tau_J = Indicator.tau_J + face.ds' * (face.penalty_geom.max(iedg) * (Data.theta*(uh_self - uh_neigh) + (1-Data.theta)*(uold_self - uold_neigh)).^2);
        Indicator.tau_N = Indicator.tau_N + face.ds' * (Data.h * (Data.theta*(grad_uh_n_self - grad_uh_n_neigh) + (1-Data.theta)*(grad_uold_n_self - grad_uold_n_neigh)).^2);
        Indicator.tau_T = Indicator.tau_T + face.ds' * (Data.h * (Data.theta*(grad_uh_t_self - grad_uh_t_neigh) + (1-Data.theta)*(grad_uold_t_self - grad_uold_t_neigh)).^2);

    end

end
