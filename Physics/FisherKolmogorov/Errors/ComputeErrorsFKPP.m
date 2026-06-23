%> @file ComputeErrorsFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Compute errors for convergence analysis
%>
%==========================================================================
%> @section classComputeErrorsFKPP description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data        Struct with problem's data
%> @param mesh	      Struct containing mesh information (region+neighbor)	
%> @param femregion   Struct containing all the information 
%>                    about the finite element approximation
%> @param c           Problem's solution
%> @param time        Current time
%>
%> @retval Error      Structure with computed L2 and dG errors 
%>                   
%==========================================================================

function [Error] = ComputeErrorsFKPP(Data, mesh, femregion, c, time)
    
    Funcs.Preallocation    = @ErrorPreallocationFKPP;
    Funcs.VolumeAssemblyST = @VolumeErrorAssemblyFKPP;
    Funcs.FacesAssembly    = @FacesErrorAssemblyFKPP; 
    Funcs.FinalMatrices    = @ErrorFKPP;

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
    
    AssembInfo.c = c;
    AssembInfo.t = time;

    [Error] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Error] = ErrorPreallocationFKPP(GenMatrices)
    Error.Volume.err_L2_loc   = GenMatrices.CellVector;
    Error.Volume.err_H1_loc   = GenMatrices.CellVector;

    Error.Faces.err_jumps_loc = GenMatrices.CellVector;
end


%% Error assembly functions
function [Error] = VolumeErrorAssemblyFKPP(Data, Error, elem, ie, id, ~, AssembInfo)
                
        % Evaluation of physical parameters
        D_ext = Data.D_ext{id}(elem.xq,elem.yq);

        % Exact solution evaluation    
        cex_loc = Data.c_ex{1}(elem.xq,elem.yq, AssembInfo.t);
        gradx_cex_loc = Data.dc_dx_ex{1}(elem.xq,elem.yq, AssembInfo.t);
        grady_cex_loc = Data.dc_dy_ex{1}(elem.xq,elem.yq, AssembInfo.t);

        % Numerical solution evaluation
        ch_loc = elem.phiq * AssembInfo.c(elem.index);
        gradx_ch_loc = elem.gradqx * AssembInfo.c(elem.index);
        grady_ch_loc = elem.gradqy * AssembInfo.c(elem.index);
        
        % Local error integral assembly
        Error.err_L2_loc = (elem.dx .* (ch_loc - cex_loc))' * (ch_loc - cex_loc); 
        Error.err_H1_loc = (elem.dx .* D_ext .* (gradx_ch_loc - gradx_cex_loc))' * (gradx_ch_loc - gradx_cex_loc) ...
                             + (elem.dx .* D_ext .* (grady_ch_loc - grady_cex_loc))' * (grady_ch_loc - grady_cex_loc);   

end

function [Error] = FacesErrorAssemblyFKPP(Data, femregion, Error, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;
    
    index_self  = sum(femregion.nbases(1:face.ie-1))+1:sum(femregion.nbases(1:face.ie));
    index_neigh = sum(femregion.nbases(1:face.neigh_ie-1))+1:sum(femregion.nbases(1:face.neigh_ie));

    % Evaluation solution and gradients
    ch_self  = face.phiedgeq * AssembInfo.c(index_self);
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        gD  = Data.DirBC{1}(face.xq,face.yq, AssembInfo.t);
        Error.err_jumps_loc = Error.err_jumps_loc + face.ds' * (face.penalty_geom.max(iedg) * (ch_self-gD).^2);

    %% Internal faces
    elseif face.neigh_ie > ie

        ch_neigh = face.phiedgeqneigh * AssembInfo.c(index_neigh);
        Error.err_jumps_loc = Error.err_jumps_loc + face.ds' * (face.penalty_geom.max(iedg) * (ch_self-ch_neigh).^2);

    end

end

%% Final error construction function
function [Error] = ErrorFKPP(Error_loc)
    Error.L2 = sqrt(sum(Error_loc.Volume.err_L2_loc));
    Error.dG = sqrt(sum(Error_loc.Volume.err_H1_loc+Error_loc.Faces.err_jumps_loc));

    Error.Cells_L2 = sqrt(Error_loc.Volume.err_L2_loc);
    Error.Cells_dG = sqrt(Error_loc.Volume.err_H1_loc+Error_loc.Faces.err_jumps_loc);
end