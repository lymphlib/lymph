%> @file  ForcingTermAssemblyPoroElastoAcoustics.m
%> @author Mattia Corti
%> @date 30 May 2026
%> @brief Choice of the assembly method for the forcing term.
%>
%==========================================================================
%> @section classForcingTermAssemblyPoroElastoAcoustics Class description
%==========================================================================
%> @brief            Choice of the assembly method for the forcing.
%>
%> @param Data          Struct with problem's data
%> @param mesh          mesh struct (region+neighbor)
%> @param femregion     Finite Element struct (see CreateDOF.m)
%> @param time          Time associated with the evaluation
%>
%> @retval F         Computed forcing term
%>                   
%==========================================================================

function [F] = ForcingTermAssemblyPoroElastoAcoustics(Data, mesh, femregion, time)
    
    Funcs.Preallocation    = @ForcingPreallocationPoroElastoAcoustics;
    Funcs.VolumeAssemblyST = @VolumeForcingAssemblyPoroElastoAcousticsST;
    Funcs.FacesAssembly    = @FacesForcingAssemblyPoroElastoAcoustics;
    Funcs.FinalMatrices    = @ForcingPoroElastoAcoustics;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = false;
    AssembInfo.assemblytrilinearforms  = false;
    
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = true;

    AssembInfo.t = time;

    AssembInfo.label = femregion.label;
    AssembInfo.MatTag = struct('Poro', 'P', ...
                                'Ela', 'E', ...
                                'Acu', 'A');

    % Control adaptivity and decide in case where to assemble 
    AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
    AssembInfo.ass_face_vec = ones(femregion.nel,1);

    % Assembly
    [F] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);

end