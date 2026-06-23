%> @file ComputeModalSolutionHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 February 2026
%> @brief Compute modal coefficient of the exact solution
%>
%==========================================================================
%> @section classComputeModalSolutionHeat description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data       Struct with problem's data
%> @param mesh       mesh struct (region+neighbor)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param time       Current time
%>
%> @retval u_mod     Modal coefficients of the exact solution 
%>                   
%==========================================================================

function [u_mod] = ComputeModalSolutionHeat(Data, mesh, femregion, time)
    
    Funcs.Preallocation    = @SolutionPreallocationHeat;
    Funcs.VolumeAssemblyST = @SolutionAssemblyHeat;
    Funcs.FinalMatrices    = @SolutionHeat;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = false;
    AssembInfo.assemblyinternalfaces   = false;
    AssembInfo.assemblytrilinearforms  = false;

    AssembInfo.computegradients        = false;
    AssembInfo.computelaplacian        = false;
    AssembInfo.computefacegradients    = false;

    AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
    AssembInfo.ass_face_vec = ones(femregion.nel,1);

    AssembInfo.t = time;
    
    [u_mod] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Solution] = SolutionPreallocationHeat(GenMatrices)
    Solution.Volume.u_loc   = GenMatrices.Vector;
end


%% Solution assembly function
function [Solution] = SolutionAssemblyHeat(Data, Solution, elem, ie, ~, nbases, AssembInfo)
        
        % Exact solution evaluation    
        uex_loc = Data.u_ex{1}(elem.xq,elem.yq,AssembInfo.t);

        % Local solution integral assembly
        Solution.u_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * uex_loc;     
end


%% Final vector construction function
function [u_mod] = SolutionHeat(Solution)
    u_mod = Solution.Volume.u_loc;
end
