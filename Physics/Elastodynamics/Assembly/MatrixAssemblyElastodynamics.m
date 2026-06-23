%> @file MatrixAssemblyElastodynamics.m
%> @author Mattia Corti
%> @date 8 May 2026
%> @brief Choice of the assembly method for the matrices.
%>
%==========================================================================
%> @section classMatrixAssemblyElastodynamics Class description
%==========================================================================
%> @brief            Choice of the assembly method for the matrices.
%>
%> @param Data                  Struct with problem's data
%> @param mesh                  mesh struct (region+neighbor)
%> @param femregion             Finite Element struct (see CreateDOF.m)
%>
%> @retval Matrices  Computed matrices
%>                   
%==========================================================================

function [Matrices] = MatrixAssemblyElastodynamics(Data, mesh, femregion)
    
    Funcs.Preallocation    = @MatrixPreallocationElastodynamics;
    Funcs.VolumeAssemblyQF = @VolumeMatricesAssemblyElastodynamicsQF;
    Funcs.VolumeAssemblyST = @VolumeMatricesAssemblyElastodynamicsST;
    Funcs.FacesAssembly    = @FacesMatricesAssemblyElastodynamics;
    Funcs.FinalMatrices    = @SIPMatricesElastodynamics;

    AssembInfo.quadrature              = Data.quadrature;
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;
    
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = true;

    AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
    AssembInfo.ass_face_vec = ones(femregion.nel,1);

    % Assembly
    [Matrices] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
