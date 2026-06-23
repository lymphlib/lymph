%> @file MatrixAssemblyHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 February 2026
%> @brief Choice of the assembly method for the matrices.
%>
%==========================================================================
%> @section classMatrixAssemblyHeat Class description
%==========================================================================
%> @brief            Choice of the assembly method for the matrices.
%>
%> @param Data                  Struct with problem's data
%> @param mesh                  mesh struct (region+neighbor)
%> @param femregion             Finite Element struct (see CreateDOF.m)
%> @param Matrices_adapt_old    Matrices of previous assembly
%>
%> @retval Matrices  Computed matrices
%>                   
%==========================================================================

function [Matrices] = MatrixAssemblyHeat(Data, mesh, femregion, Matrices_adapt_old)
    
    Funcs.Preallocation    = @MatrixPreallocationHeat;
    Funcs.VolumeAssemblyQF = @VolumeMatricesAssemblyHeatQF;
    Funcs.VolumeAssemblyST = @VolumeMatricesAssemblyHeatST;
    Funcs.FacesAssembly    = @FacesMatricesAssemblyHeat;
    Funcs.FinalMatrices    = @SIPMatricesHeat;

    AssembInfo.quadrature              = Data.quadrature;
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;
    
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = true;

    % Control adaptivity and decide in case where to assemble 
    if not(isempty(Matrices_adapt_old))
        AssembInfo.Matrices_adapt_old = Matrices_adapt_old.Adaptivity;
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
    [Matrices] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
