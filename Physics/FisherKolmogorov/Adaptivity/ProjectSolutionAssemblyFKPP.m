%> @file   ProjectSolutionAssemblyFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 15 April 2026
%> @brief Projection of the solution from the old femregion to the new one.
%>
%==========================================================================
%> @section classProjectionSolutionAssemblyFKPP Class description
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

function [SolutionAdapt] = ProjectSolutionAssemblyFKPP(Data, mesh, femregion, time, Solution)
    
    Funcs.Preallocation    = @ProjectionPreallocation;
    Funcs.VolumeAssemblyST = @VolumeProjectionAssemblyST;
    Funcs.FinalMatrices    = @ForcingFKPP;

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
    AssembInfo.Matrices_adapt_old.c_loc     = Solution.c_h;
    AssembInfo.Matrices_adapt_old.c_old_loc = Solution.c_old;
    if isfield(Solution,"c_oold")
        AssembInfo.Matrices_adapt_old.c_oold_loc = Solution.c_oold;
    end

    AssembInfo.ass_vol_vec  = ones(length(femregion.degree),1);

    % Assembly
    [SolutionAdapt] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
    SolutionAdapt.Indicator = Solution.Indicator;

    if not(isfield(Solution,"c_oold"))
        SolutionAdapt = rmfield(SolutionAdapt,"c_oold");
    end
end


%% Preallocation of the projection vector
function [Solution] = ProjectionPreallocation(GenMatrices)
    Solution.Volume.c_loc      = GenMatrices.Vector;
    Solution.Volume.c_old_loc  = GenMatrices.Vector;
    Solution.Volume.c_oold_loc = GenMatrices.Vector;
end

%% Assembly of the projection
function [Solution] = VolumeProjectionAssemblyST(Data, Solution, elem, ie, id, nbases, AssembInfo)
        
        c_h_loc   = elem.phiq_old*AssembInfo.Matrices_adapt_old.c_loc(elem.index_old);
        c_old_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.c_old_loc(elem.index_old);

        c_h_loc   = (elem.dx .* elem.phiq)' * c_h_loc;
        c_old_loc = (elem.dx .* elem.phiq)' * c_old_loc;

        M_loc = (elem.dx .* elem.phiq)' * elem.phiq;

        Solution.c_loc(1:nbases,1)     = M_loc\c_h_loc;
        Solution.c_old_loc(1:nbases,1) = M_loc\c_old_loc;

        if isfield(AssembInfo.Matrices_adapt_old,"c_oold_loc")
            c_oold_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.c_oold_loc(elem.index_old);
            c_oold_loc = (elem.dx .* elem.phiq)' * c_oold_loc;
            Solution.c_oold_loc(1:nbases,1) = M_loc\c_oold_loc;
        end

end

%% Final solution assembly
function [Solution] = ForcingFKPP(Solution_loc)
    Solution.c_h     = Solution_loc.Volume.c_loc;
    Solution.c_old   = Solution_loc.Volume.c_old_loc;
    Solution.c_oold  = Solution_loc.Volume.c_oold_loc;
end

