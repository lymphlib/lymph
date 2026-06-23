%> @file ComputeModalSolutionFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Compute modal coefficient of the exact solution
%>
%==========================================================================
%> @section classComputeModalSolutionFHN description
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

function [Solution] = ComputeModalSolutionFHN(Data, mesh, femregion, time)
    
    Funcs.Preallocation    = @SolutionPreallocationFHN;
    Funcs.VolumeAssemblyST = @SolutionAssemblyFHN;
    Funcs.FinalMatrices    = @SolutionFHN;

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
    
    [Solution] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Solution] = SolutionPreallocationFHN(GenMatrices)
    Solution.Volume.u_loc   = GenMatrices.Vector;
    Solution.Volume.w_loc   = GenMatrices.Vector;
end


%% Solution assembly function
function [Solution] = SolutionAssemblyFHN(Data, Solution, elem, ie, ~, nbases, AssembInfo)
        
        % Exact solution evaluation    
        uex_loc = Data.u_ex{1}(elem.xq,elem.yq,AssembInfo.t);
        wex_loc = Data.w_ex{1}(elem.xq,elem.yq,AssembInfo.t);

        % Local solution integral assembly
        Solution.u_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * uex_loc;     
        Solution.w_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * wex_loc;     
end


%% Final vector construction function
function [SolutionFinal] = SolutionFHN(Solution)
    SolutionFinal.u_h = Solution.Volume.u_loc;
    SolutionFinal.w_h = Solution.Volume.w_loc;
end
