%> @file  IPDGMatrixAssemblyFKPP.m
%> @author Mattia Corti
%> @date 15 September 2025
%> @brief Choice of the assembly properties for the matrices of IP method.
%>
%==========================================================================
%> @section classIPDGMatrixAssemblyFKPP Class description
%==========================================================================
%> @brief            Choice of the assembly method for the matrices.
%>
%> @param Data       		Struct with problem's data
%> @param mesh       		Mesh struct (region+neighbor)
%> @param femregion  		Finite Element struct (see CreateDOF.m)
%> @param Matrices_adapt_old 	Matrices associated with the previous adaptivity iteration
%>
%> @retval Matrices  Computed matrices
%>                   
%==========================================================================

function [Matrices] = IPDGMatrixAssemblyFKPP(Data, mesh, femregion, Matrices_adapt_old)
    
    Funcs.Preallocation      = @IPMatrixPreallocationFK;
    Funcs.VolumeAssemblyQF   = @IPVolumeMatricesAssemblyFKQF;
    Funcs.VolumeAssemblyST   = @IPVolumeMatricesAssemblyFKST;
    Funcs.Volume3LAssemblyQF = @IPVolume3LMatricesAssemblyFKQF;
    Funcs.Volume3LAssemblyST = @IPVolume3LMatricesAssemblyFKST;
    Funcs.FacesAssembly      = @IPFacesMatricesAssemblyFK;
    Funcs.FinalMatrices      = @SIPMatricesFK;

    AssembInfo.quadrature              = Data.quadrature;
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = true;

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

    if ~Data.isotropy
        AssembInfo.Tens = readtable(fullfile(Data.foldername,Data.AxnDiffFile));
    end

    [Matrices] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
