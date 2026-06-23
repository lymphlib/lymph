%> @file MatrixAssemblyPoroElastoAcoustics.m
%> @author Mattia Corti
%> @date 29 May 2026
%> @brief Choice of the assembly method for the matrices.
%>
%==========================================================================
%> @section classMatrixAssemblyPoroElastoAcoustics Class description
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

function [Matrices] = MatrixAssemblyPoroElastoAcoustics(Data, mesh, femregion)
    
    Funcs.Preallocation    = @MatrixPreallocationPoroElastoAcoustics;
    Funcs.VolumeAssemblyQF = @VolumeMatricesAssemblyPoroElastoAcousticsQF;
    Funcs.VolumeAssemblyST = @VolumeMatricesAssemblyPoroElastoAcousticsST;
    Funcs.FacesAssembly    = @FacesMatricesAssemblyPoroElastoAcoustics;
    Funcs.FinalMatrices    = @SIPMatricesPoroElastoAcoustics;

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

    AssembInfo.label = femregion.label;
    AssembInfo.MatTag = struct('Poro',   'PP','Ela',   'EE','Acu',    'AA',...
                               'PoroAcu','PA','ElaAcu','EA','PoroEla','PE');

    % Assembly
    [Matrices] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
