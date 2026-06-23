%> @file   ProjectSolutionAssemblyHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 18 February 2026
%> @brief Projection of the solution from the old femregion to the new one.
%>
%==========================================================================
%> @section classProjectionSolutionAssemblyHeat Class description
%==========================================================================
%> @brief            Projection of the solution from the old femregion to the new one.
%>
%> @param Data          Struct with problem's data
%> @param mesh          mesh struct (region+neighbor)
%> @param femregion     Finite Element struct (see CreateDOF.m)
%> @param time          Time associated with the evaluation
%> @param Solution      Solution to project
%>
%> @retval SolutionAdapt Solution projected on the new femregion
%>                   
%==========================================================================

function [SolutionAdapt] = ProjectSolutionAssemblyHeat(Data, mesh, femregion, time, Solution)
    
    Funcs.Preallocation    = @ProjectionPreallocation;
    Funcs.VolumeAssemblyST = @VolumeProjectionAssemblyST;
    Funcs.FinalMatrices    = @ForcingHeat;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = false;
    AssembInfo.assemblyinternalfaces   = false;
    AssembInfo.assemblytrilinearforms  = false;
    
    AssembInfo.computegradients        = false;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = true;
    AssembInfo.computeoldbases         = true;

    AssembInfo.t = time;

    % Control adaptivity and decide in case where to assemble
    AssembInfo.Matrices_adapt_old.u_loc     = Solution.u_h;
    AssembInfo.Matrices_adapt_old.u_old_loc = Solution.u_old;

    AssembInfo.ass_vol_vec  = ones(length(femregion.degree),1);

    % Assembly
    [SolutionAdapt] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
    SolutionAdapt.Indicator = Solution.Indicator;
end


%% Preallocation of the projection vector
function [Solution] = ProjectionPreallocation(GenMatrices)
    Solution.Volume.u_loc      = GenMatrices.Vector;
    Solution.Volume.u_old_loc  = GenMatrices.Vector;
end

%% Assembly of the projection
function [Solution] = VolumeProjectionAssemblyST(Data, Solution, elem, ie, id, nbases, AssembInfo)
        
        u_h_loc   = elem.phiq_old*AssembInfo.Matrices_adapt_old.u_loc(elem.index_old);
        u_old_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.u_old_loc(elem.index_old);

        U_h_loc   = (elem.dx .* elem.phiq)' * u_h_loc;
        U_old_loc = (elem.dx .* elem.phiq)' * u_old_loc;

        M_loc = (elem.dx .* elem.phiq)' * elem.phiq;

        Solution.u_loc(1:nbases,1)     = M_loc\U_h_loc;
        Solution.u_old_loc(1:nbases,1) = M_loc\U_old_loc;
end

%% Final solution assembly
function [Solution] = ForcingHeat(Solution_loc)
    Solution.u_h     = Solution_loc.Volume.u_loc;
    Solution.u_old   = Solution_loc.Volume.u_old_loc;
end

