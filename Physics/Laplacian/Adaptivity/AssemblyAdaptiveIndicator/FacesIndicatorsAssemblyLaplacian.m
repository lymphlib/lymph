%> @file   FacesIndicatorsAssemblyLaplacian.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   12 May 2026
%> @brief Computation of the face indicators for Laplacian problem.
%>
%==========================================================================
%> @section classFacesIndicatorsAssemblyLaplacian Class description
%==========================================================================
%> @brief            Computation of the face indicators for Laplacian problem.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Indicator  Indicator struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Indicator  Indicator struct
%>                   
%==========================================================================

function [Indicator] = FacesIndicatorsAssemblyLaplacian(Data, femregion, Indicator, face, AssembInfo)
   
    iedg = face.iedg;
    
    index_self  = sum(femregion.nbases(1:face.ie-1))+1:sum(femregion.nbases(1:face.ie));
    index_neigh = sum(femregion.nbases(1:face.neigh_ie-1))+1:sum(femregion.nbases(1:face.neigh_ie));

    % Evaluation solution and gradients
    uh_self  = face.phiedgeq * AssembInfo.uh(index_self);
    grad_uh_n_self  =  (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy) * AssembInfo.uh(index_self);
    grad_uh_t_self  = (-face.ny * face.gradedgeqx + face.nx * face.gradedgeqy) * AssembInfo.uh(index_self);
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        gD  = Data.DirBC{1}(face.xq,face.yq);
        gDx = Data.gradDirBC{1}(face.xq,face.yq);
        gDy = Data.gradDirBC{2}(face.xq,face.yq);

        Indicator.tau_J = Indicator.tau_J + face.ds' * (face.penalty_geom.max(iedg) * (uh_self-gD).^2);
        Indicator.tau_T = Indicator.tau_T + face.ds' * (Data.h * (grad_uh_t_self - (-face.ny*gDx + face.nx*gDy)).^2);
   
    %% Internal faces
    elseif face.neigh_ie > 0

        uh_neigh        = face.phiedgeqneigh * AssembInfo.uh(index_neigh);
        grad_uh_n_neigh =  (face.nx * face.gradedgeqxneigh + face.ny * face.gradedgeqyneigh) * AssembInfo.uh(index_neigh);
        grad_uh_t_neigh = (-face.ny * face.gradedgeqxneigh + face.nx * face.gradedgeqyneigh) * AssembInfo.uh(index_neigh);
    
        Indicator.tau_J = Indicator.tau_J + face.ds' * (face.penalty_geom.max(iedg) * (uh_self-uh_neigh).^2);
        Indicator.tau_N = Indicator.tau_N + face.ds' * (Data.h * (grad_uh_n_self - grad_uh_n_neigh).^2);
        Indicator.tau_T = Indicator.tau_T + face.ds' * (Data.h * (grad_uh_t_self - grad_uh_t_neigh).^2);

    end

end