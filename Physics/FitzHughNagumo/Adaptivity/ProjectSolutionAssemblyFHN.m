%> @file   ProjectSolutionAssemblyFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Projection of the solution from the old femregion to the new one.
%>
%==========================================================================
%> @section classProjectionSolutionAssemblyFHN Class description
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

function [SolutionAdapt] = ProjectSolutionAssemblyFHN(Data, mesh, femregion, time, Solution)
    
    Funcs.Preallocation    = @ProjectionPreallocation;
    Funcs.VolumeAssemblyST = @VolumeProjectionAssemblyST;
    Funcs.FinalMatrices    = @ForcingFHN;

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
    AssembInfo.Matrices_adapt_old.u_loc     = Solution.t.u_h;
    AssembInfo.Matrices_adapt_old.u_old_loc = Solution.t_old.u_h;

    AssembInfo.Matrices_adapt_old.w_loc     = Solution.t.w_h;
    AssembInfo.Matrices_adapt_old.w_old_loc = Solution.t_old.w_h;

    if isfield(Solution,"t_oold")
        AssembInfo.Matrices_adapt_old.u_oold_loc = Solution.t_oold.u_h;
        AssembInfo.Matrices_adapt_old.w_oold_loc = Solution.t_oold.w_h;
    end

    AssembInfo.ass_vol_vec  = ones(length(femregion.degree),1);

    % Assembly
    [SolutionAdapt] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
    
    SolutionAdapt.t.Indicator = Solution.t.Indicator;

end


%% Preallocation of the projection vector
function [Solution] = ProjectionPreallocation(GenMatrices)
    Solution.Volume.u_loc      = GenMatrices.Vector;
    Solution.Volume.u_old_loc  = GenMatrices.Vector;
    Solution.Volume.u_oold_loc = GenMatrices.Vector;

    Solution.Volume.w_loc      = GenMatrices.Vector;
    Solution.Volume.w_old_loc  = GenMatrices.Vector;
    Solution.Volume.w_oold_loc = GenMatrices.Vector;
end

%% Assembly of the projection
function [Solution] = VolumeProjectionAssemblyST(Data, Solution, elem, ie, id, nbases, AssembInfo)
        
        u_h_loc   = elem.phiq_old*AssembInfo.Matrices_adapt_old.u_loc(elem.index_old);
        u_old_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.u_old_loc(elem.index_old);
        u_h_loc   = (elem.dx .* elem.phiq)' * u_h_loc;
        u_old_loc = (elem.dx .* elem.phiq)' * u_old_loc;
        
        w_h_loc   = elem.phiq_old*AssembInfo.Matrices_adapt_old.w_loc(elem.index_old);
        w_old_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.w_old_loc(elem.index_old);
        w_h_loc   = (elem.dx .* elem.phiq)' * w_h_loc;
        w_old_loc = (elem.dx .* elem.phiq)' * w_old_loc;

        M_loc = (elem.dx .* elem.phiq)' * elem.phiq;

        Solution.u_loc(1:nbases,1)     = M_loc\u_h_loc;
        Solution.u_old_loc(1:nbases,1) = M_loc\u_old_loc;
        
        Solution.w_loc(1:nbases,1)     = M_loc\w_h_loc;
        Solution.w_old_loc(1:nbases,1) = M_loc\w_old_loc;
        
        if isfield(AssembInfo.Matrices_adapt_old,'u_oold_loc')

            u_oold_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.u_oold_loc(elem.index_old);
            u_oold_loc = (elem.dx .* elem.phiq)' * u_oold_loc;
            w_oold_loc = elem.phiq_old*AssembInfo.Matrices_adapt_old.w_oold_loc(elem.index_old);
            w_oold_loc = (elem.dx .* elem.phiq)' * w_oold_loc;

            Solution.u_oold_loc(1:nbases,1) = M_loc\u_oold_loc;
            Solution.w_oold_loc(1:nbases,1) = M_loc\w_oold_loc;

        end
        
end

%% Final solution assembly
function [Solution] = ForcingFHN(Solution_loc)
    Solution.t.u_h      = Solution_loc.Volume.u_loc;
    Solution.t_old.u_h  = Solution_loc.Volume.u_old_loc;
    Solution.t_oold.u_h = Solution_loc.Volume.u_oold_loc;

    Solution.t.w_h      = Solution_loc.Volume.w_loc;
    Solution.t_old.w_h  = Solution_loc.Volume.w_old_loc;
    Solution.t_oold.w_h = Solution_loc.Volume.w_oold_loc;
end

