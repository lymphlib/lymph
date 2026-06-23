%> @file  ForcingTermAssemblyHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 February 2026
%> @brief Choice of the assembly method for the forcing term.
%>
%==========================================================================
%> @section classForcingTermAssemblyHeat Class description
%==========================================================================
%> @brief            Choice of the assembly method for the forcing.
%>
%> @param Data          Struct with problem's data
%> @param mesh          mesh struct (region+neighbor)
%> @param femregion     Finite Element struct (see CreateDOF.m)
%> @param time          Time associated with the evaluation
%> @param F_adapt_old   Forcing term of previous assembly
%>
%> @retval F         Computed forcing term
%>                   
%==========================================================================

function [F] = ForcingTermAssemblyHeat(Data, mesh, femregion, time, F_adapt_old)
    
    Funcs.Preallocation    = @ForcingPreallocationHeat;
    Funcs.VolumeAssemblyST = @VolumeForcingAssemblyHeatST;
    Funcs.FacesAssembly    = @FacesForcingAssemblyHeat;
    Funcs.FinalMatrices    = @ForcingHeat;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = (Data.TagApplyBCs == 1);
    AssembInfo.assemblyinternalfaces   = false;
    AssembInfo.assemblytrilinearforms  = false;
    
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = true;

    AssembInfo.t = time;

    % Control adaptivity and decide in case where to assemble 
    if not(isempty(F_adapt_old))
        AssembInfo.Matrices_adapt_old = F_adapt_old.Adaptivity;
        AssembInfo.ass_vol_vec  = (femregion.degree>femregion.degree_old);
        AssembInfo.ass_face_vec = (femregion.degree~=femregion.degree_old);
        neighs_upd = unique([mesh.neighbor.neigh{AssembInfo.ass_face_vec}]);
        neighs_upd(neighs_upd<=0) = [];
        AssembInfo.ass_face_vec(neighs_upd) = 1;
    else
        AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
        AssembInfo.ass_face_vec = ones(femregion.nel,1);
    end

    % Assembly
    [F] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);

end