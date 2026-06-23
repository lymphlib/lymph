%> @file  ForcingTermAssemblyFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Choice of the assembly method for the forcing term.
%>
%==========================================================================
%> @section classForcingTermAssemblyFHN Class description
%==========================================================================
%> @brief            Choice of the assembly method for the matrices.
%>
%> @param Data        Struct with problem's data
%> @param mesh        mesh struct (region+neighbor)
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param time        Time for the evaluation of the term
%> @param Solution    Solution structure
%> @param F_adapt_old Forcing term associated with the previous adaptivity step
%>
%> @retval F          Computed forcing term
%>
%==========================================================================

function [F] = ForcingTermAssemblyFHN(Data, mesh, femregion, time, Solution, F_adapt_old)

Funcs.Preallocation     = @ForcingPreallocationFHN;
Funcs.VolumeAssemblyST  = @VolumeForcingAssemblyFHNST;
Funcs.FacesAssembly     = @FacesForcingAssemblyFHN;
Funcs.FinalMatrices     = @ForcingFHN;

AssembInfo.quadrature              = "ST";
AssembInfo.assemblyvolume          = true;
AssembInfo.assemblyfaces           = (Data.TagApplyBCs == 1);
AssembInfo.assemblyinternalfaces   = false;
AssembInfo.assemblytrilinearforms  = false;

AssembInfo.computegradients        = true;
AssembInfo.computelaplacian        = false;
AssembInfo.computefacegradients    = true;

AssembInfo.t                       = time;

AssembInfo.u                       = Solution.u_h;
AssembInfo.w                       = Solution.w_h;

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

end