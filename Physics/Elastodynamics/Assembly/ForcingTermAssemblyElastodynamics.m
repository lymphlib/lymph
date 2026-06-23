%> @file  ForcingTermAssemblyElastodynamics.m
%> @author Mattia Corti
%> @date 8 May 2026
%> @brief Choice of the assembly method for the forcing term.
%>
%==========================================================================
%> @section classForcingTermAssemblyElastodynamics Class description
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

function [F] = ForcingTermAssemblyElastodynamics(Data, mesh, femregion, time)
    
    Funcs.Preallocation    = @ForcingPreallocationElastodynamics;
    Funcs.VolumeAssemblyST = @VolumeForcingAssemblyElastodynamicsST;
    Funcs.FacesAssembly    = @FacesForcingAssemblyElastodynamics;
    Funcs.FinalMatrices    = @ForcingElastodynamics;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = false;
    AssembInfo.assemblytrilinearforms  = false;
    
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = true;

    AssembInfo.t = time;

    % Control adaptivity and decide in case where to assemble 
    AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
    AssembInfo.ass_face_vec = ones(femregion.nel,1);

    % Assembly
    [F] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);

end