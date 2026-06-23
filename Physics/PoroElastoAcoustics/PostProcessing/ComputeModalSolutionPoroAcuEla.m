%> @file ComputeModalSolutionPoroAcuEla.m
%> @author Mattia Corti
%> @date 5 June 2026
%> @brief Compute modal coefficient of the exact solution
%>
%==========================================================================
%> @section classComputeModalSolutionPoroAcuEla description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data       Struct with problem's data
%> @param mesh       mesh struct (region+neighbor)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param time       Time at which the solution is evaluated
%>
%> @retval Solution  Structure with modal coefficients of the exact solution
%>                   
%==========================================================================

function [Solution] = ComputeModalSolutionPoroAcuEla(Data, mesh, femregion, time)
    
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
    
    AssembInfo.label = femregion.label;
    AssembInfo.MatTag = struct('Poro', 'P', ...
                                'Ela', 'E', ...
                                'Acu', 'A');

    [Solution] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end


%% Preallocation function
function [Solution] = SolutionPreallocationElastodynamics(GenMatrices)

    % Preallocation poroelasticity solution
    Solution.Volume.Poro.up1_loc  = GenMatrices.Vector;
    Solution.Volume.Poro.up2_loc  = GenMatrices.Vector;
    Solution.Volume.Poro.wp1_loc  = GenMatrices.Vector;
    Solution.Volume.Poro.wp2_loc  = GenMatrices.Vector;

    % Preallocation elasticity solution
    Solution.Volume.Ela.ue1_loc   = GenMatrices.Vector;
    Solution.Volume.Ela.ue2_loc   = GenMatrices.Vector;

    % Preallocation acoustics solution
    Solution.Volume.Acu.phi_a_loc = GenMatrices.Vector;

end


%% Solution assembly function
function [Solution] = SolutionAssemblyElastodynamics(Data, Solution, elem, ie, ~, nbases, AssembInfo)
        
        switch AssembInfo.label(ie)
            case 'P'
                % Exact solution evaluation                
                up1ex_loc = Data.up_ex{1}(elem.xq,elem.yq);
                up2ex_loc = Data.up_ex{2}(elem.xq,elem.yq);
                wp1ex_loc = Data.wp_ex{1}(elem.xq,elem.yq);
                wp2ex_loc = Data.wp_ex{2}(elem.xq,elem.yq);
                    
                % Local solution integral assembly          
                Solution.Poro.up1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * up1ex_loc; 
                Solution.Poro.up2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * up2ex_loc; 
                Solution.Poro.wp1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * wp1ex_loc; 
                Solution.Poro.wp2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * wp2ex_loc;     
                            
            case 'E'
                % Exact solution evaluation                
                ue1ex_loc = Data.ue_ex{1}(elem.xq,elem.yq);
                ue2ex_loc = Data.ue_ex{2}(elem.xq,elem.yq);
                    
                % Local solution integral assembly          
                Solution.Ela.ue1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * ue1ex_loc; 
                Solution.Ela.ue2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * ue2ex_loc; 

            case 'A'
                % Exact solution evaluation                
                phiex_loc = Data.phi_ex{1}(elem.xq,elem.yq);
                
                % Local solution integral assembly          
                Solution.Acu.phi_a_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * phiex_loc; 
        end
end


%% Final vector construction function
function [Solution] = SolutionElastodynamics(Solution_loc)
    
    % Poroelasticity
    Solution.up = [Solution_loc.Volume.Poro.up1_loc; Solution_loc.Volume.Poro.up2_loc];
    Solution.wp = [Solution_loc.Volume.Poro.wp1_loc; Solution_loc.Volume.Poro.wp2_loc];
    
    % Elastodynamics
    Solution.ue = [Solution_loc.Volume.Ela.ue1_loc; Solution_loc.Volume.Ela.ue2_loc];

    % Acoustics
    Solution.phi = Solution_loc.Volume.Acu.phi_a_loc;
    
end
