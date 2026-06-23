%> @file ComputeErrorsHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 3 February 2026
%> @brief Compute errors for convergence analysis
%>
%==========================================================================
%> @section classComputeErrorsHeat description
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

function [Error] = ComputeErrorsHeat(Data, mesh, femregion, U, time)
    
    Funcs.Preallocation    = @ErrorPreallocationLaplacian;
    Funcs.VolumeAssemblyST = @VolumeErrorAssemblyLaplacian;
    Funcs.FacesAssembly    = @FacesErrorAssemblyLaplacian; 
    Funcs.FinalMatrices    = @ErrorLaplacian;

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
    
    AssembInfo.U = U;
    AssembInfo.t = time;

    [Error] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
    Error.NDoF = femregion.ndof;
   
end


%% Preallocation function
function [Error] = ErrorPreallocationLaplacian(GenMatrices)
    Error.Volume.err_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_H1_loc   = GenMatrices.CellVector;

    Error.Faces.err_jumps_loc = GenMatrices.CellVector;
end


%% Error assembly functions
function [Error] = VolumeErrorAssemblyLaplacian(Data, Error, elem, ie, id, ~, AssembInfo)
                
        % Evaluation of physical parameters
        mu = Data.mu{id}(elem.xq,elem.yq);

        % Exact solution evaluation    
        uex_loc = Data.u_ex{1}(elem.xq,elem.yq, AssembInfo.t);
        gradx_uex_loc = Data.du_dx_ex{1}(elem.xq,elem.yq, AssembInfo.t);
        grady_uex_loc = Data.du_dy_ex{1}(elem.xq,elem.yq, AssembInfo.t);

        % Numerical solution evaluation
        uh_loc = elem.phiq * AssembInfo.U(elem.index);
        gradx_uh_loc = elem.gradqx * AssembInfo.U(elem.index);
        grady_uh_loc = elem.gradqy * AssembInfo.U(elem.index);
        
        % Local error integral assembly
        Error.err_L2_loc = (elem.dx .* (uh_loc - uex_loc))' * (uh_loc - uex_loc); 
        Error.err_H1_loc = (elem.dx .* mu .* (gradx_uh_loc - gradx_uex_loc))' * (gradx_uh_loc - gradx_uex_loc) ...
                             + (elem.dx .* mu .* (grady_uh_loc - grady_uex_loc))' * (grady_uh_loc - grady_uex_loc);   

end

function [Error] = FacesErrorAssemblyLaplacian(Data, femregion, Error, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;
    
    index_self  = sum(femregion.nbases(1:face.ie-1))+1:sum(femregion.nbases(1:face.ie));
    index_neigh = sum(femregion.nbases(1:face.neigh_ie-1))+1:sum(femregion.nbases(1:face.neigh_ie));

    % Evaluation solution and gradients
    uh_self  = face.phiedgeq * AssembInfo.U(index_self);
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        gD  = Data.DirBC{1}(face.xq,face.yq, AssembInfo.t);
        Error.err_jumps_loc = Error.err_jumps_loc + face.ds' * (face.penalty_geom.max(iedg) * (uh_self-gD).^2);

    %% Internal faces
    elseif face.neigh_ie > ie

        uh_neigh = face.phiedgeqneigh * AssembInfo.U(index_neigh);
        Error.err_jumps_loc = Error.err_jumps_loc + face.ds' * (face.penalty_geom.max(iedg) * (uh_self-uh_neigh).^2);

    end

end

%% Final error construction function
function [Error] = ErrorLaplacian(Error_loc)
    Error.L2 = sqrt(sum(Error_loc.Volume.err_L2_loc));
    Error.dG = sqrt(sum(Error_loc.Volume.err_H1_loc+Error_loc.Faces.err_jumps_loc));

    Error.Cells_L2 = sqrt(Error_loc.Volume.err_L2_loc);
    Error.Cells_dG = sqrt(Error_loc.Volume.err_H1_loc+Error_loc.Faces.err_jumps_loc);
end
