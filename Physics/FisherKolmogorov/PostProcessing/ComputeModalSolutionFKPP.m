%> @file ComputeModalSolutionFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 14 April 2026
%> @brief Compute modal coefficient of the exact solution
%>
%==========================================================================
%> @section classComputeModalSolutionFKPP description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data       Struct with problem's data
%> @param mesh       mesh struct (region+neighbor)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param time       Current time
%>
%> @retval c_mod     Modal coefficients of the exact solution 
%>                   
%==========================================================================

function [c_mod] = ComputeModalSolutionFKPP(Data, mesh, femregion, time)
    
    Funcs.Preallocation    = @SolutionPreallocationFKPP;
    Funcs.VolumeAssemblyST = @SolutionAssemblyFKPP;
    Funcs.FinalMatrices    = @SolutionFKPP;

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
    
    [c_mod] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Solution] = SolutionPreallocationFKPP(GenMatrices)
    Solution.Volume.c_loc   = GenMatrices.Vector;
end


%% Solution assembly function
function [Solution] = SolutionAssemblyFKPP(Data, Solution, elem, ie, ~, nbases, AssembInfo)
        
        % Exact solution evaluation    
        cex_loc = Data.c_ex{1}(elem.xq,elem.yq,AssembInfo.t);

        % Local solution integral assembly
        Solution.c_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * cex_loc;     
end


%% Final vector construction function
function [c_mod] = SolutionFKPP(Solution)
    c_mod = Solution.Volume.c_loc;
end
