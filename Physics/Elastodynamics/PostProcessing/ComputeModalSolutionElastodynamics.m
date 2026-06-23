%> @file ComputeModalSolutionElastodynamics.m
%> @author Mattia Corti
%> @date 8 May 2026
%> @brief Compute modal coefficient of the exact solution
%>
%==========================================================================
%> @section classComputeModalSolutionElastodynamics description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data       Struct with problem's data
%> @param mesh       mesh struct (region+neighbor)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param time       Time at which the solution is evaluated
%>
%> @retval u_mod     Modal coefficients of the exact solution (displacement
%                    and velocity)
%>                   
%==========================================================================

function [u_mod] = ComputeModalSolutionElastodynamics(Data, mesh, femregion, time)
    
    Funcs.Preallocation    = @SolutionPreallocationElastodynamics;
    Funcs.VolumeAssemblyST = @SolutionAssemblyElastodynamics;
    Funcs.FinalMatrices    = @SolutionElastodynamics;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = false;
    AssembInfo.assemblyinternalfaces   = false;
    AssembInfo.assemblytrilinearforms  = false;

    AssembInfo.computegradients        = false;
    AssembInfo.computeElastodynamics   = false;
    AssembInfo.computefacegradients    = false;

    AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
    AssembInfo.ass_face_vec = ones(femregion.nel,1);

    AssembInfo.t = time;
    
    [u_mod] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Solution] = SolutionPreallocationElastodynamics(GenMatrices)
    Solution.Volume.ue1_loc   = GenMatrices.Vector;
    Solution.Volume.ue2_loc   = GenMatrices.Vector;
    Solution.Volume.ve1_loc   = GenMatrices.Vector;
    Solution.Volume.ve2_loc   = GenMatrices.Vector;
end


%% Solution assembly function
function [Solution] = SolutionAssemblyElastodynamics(Data, Solution, elem, ie, ~, nbases, AssembInfo)
        
        % Exact solution evaluation                
        ue1ex_loc = Data.ue_ex{1}(elem.xq,elem.yq,AssembInfo.t);
        ue2ex_loc = Data.ue_ex{2}(elem.xq,elem.yq,AssembInfo.t);
        ve1ex_loc = Data.due_dt_ex{1}(elem.xq,elem.yq,AssembInfo.t);
        ve2ex_loc = Data.due_dt_ex{2}(elem.xq,elem.yq,AssembInfo.t);
            
        % Local solution integral assembly          
        Solution.ue1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * ue1ex_loc; 
        Solution.ue2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * ue2ex_loc; 
        Solution.ve1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * ve1ex_loc; 
        Solution.ve2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * ve2ex_loc;     
end


%% Final vector construction function
function [u_mod] = SolutionElastodynamics(Solution_loc)
    u_mod = [Solution_loc.Volume.ue1_loc;
             Solution_loc.Volume.ue2_loc;
             Solution_loc.Volume.ve1_loc;
             Solution_loc.Volume.ve2_loc];
end
