%> @file  ForcingTermAssemblyFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Choice of the assembly method for the forcing term.
%>
%==========================================================================
%> @section classForcingTermAssemblyFKPP Class description
%==========================================================================
%> @brief            Choice of the assembly method for the matrices.
%>
%> @param Data        Struct with problem's data
%> @param mesh        Mesh struct (region+neighbor)
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param time        Time for the evaluation of the term
%> @param F_adapt_old Forcing terms associated with the previous adaptivity iteration
%>
%> @retval F_new     Computed forcing term
%>                   
%==========================================================================

function [F] = ForcingTermAssemblyFKPP(Data, mesh, femregion, time, F_adapt_old)
    
    if ~Data.homog_source_f || (Data.TagApplyBCs == 1)
    
        Funcs.Preallocation     = @ForcingPreallocationFK;
        Funcs.TensorFileImport  = @TensorFileImportFK;
        Funcs.VolumeAssemblyST  = @VolumeForcingAssemblyFKST;
        Funcs.FacesAssembly     = @FacesForcingAssemblyFK;
        Funcs.FinalMatrices     = @ForcingFK;
        
        AssembInfo.quadrature              = "ST";
        AssembInfo.assemblyvolume          = ~Data.homog_source_f;
        AssembInfo.assemblyfaces           = (Data.TagApplyBCs == 1);
        AssembInfo.assemblyinternalfaces   = false;
        AssembInfo.assemblytrilinearforms  = false;
        
        AssembInfo.computegradients        = true;
        AssembInfo.computelaplacian        = false;
        AssembInfo.computefacegradients    = true;
        
        AssembInfo.t                       = time;
        
        if ~Data.isotropy
            AssembInfo.Tens = readtable(fullfile(Data.foldername,Data.AxnDiffFile));
        end
        
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
        
        [F] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
    
    else
        F.F = zeros(sum(femregion.nbases),1);
    end
end
