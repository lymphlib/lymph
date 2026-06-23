%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti, Caterina Leimer Saglio
%> @date 8 May 2026
%> @brief Compute the solution for post-processing for the elastodynamics equation.
%>
%==========================================================================
%> @section classElastodynamicsGetSolutionPostProcessing Class description
%==========================================================================
%> @brief  Compute the solution for post-processing for the heat equation.
%>
%> @param Data        Struct with problem's data
%> @param femregion   Struct for finite elements data
%> @param neighbor    Struct for neighboring elements data
%> @param Solution    Modal solution and adaptive indicators
%> @param t           Current time
%>
%> @retval Xh         Cell array containing a struct with approximate
%>                    solutions strings useful for post processing
%> @retval Xexact     Cell array containing a struct with exact
%>                    solutions strings useful for post processing
%>
%==========================================================================

function [Xh, Xexact] = GetSolutionPostProcessing(Data, femregion, neighbor, Solution, t)

    SolutionInfo.outputnames = {'u','v','p'};
    SolutionInfo.outputsizes = { 2 , 2 , 1 };
    SolutionInfo.phys_param  = {@(xq, yq, t, id) Data.lam_el{id}(xq,yq)};
    SolutionInfo.outputsol   = { [Solution(1:femregion.ndof),Solution(femregion.ndof+1:2*femregion.ndof)], ...
                                 [Solution(2*femregion.ndof+1:3*femregion.ndof),Solution(3*femregion.ndof+1:4*femregion.ndof)], ...
                                 [Solution(1:femregion.ndof),Solution(femregion.ndof+1:2*femregion.ndof)]};
    SolutionInfo.outputexpr  = {@(phi, gradx, grady, u) [phi*u(:,1), phi*u(:,2)], ...
                                @(phi, gradx, grady, v) [phi*v(:,1), phi*v(:,2)], ...
                                @(phi, gradx, grady, u, lam) -lam.*(gradx*u(:,1)+grady*u(:,2))};
    SolutionInfo.outputexact = {@(xq, yq, t) [Data.ue_ex{1}(xq, yq, t), Data.ue_ex{2}(xq, yq, t)], ...
                                @(xq, yq, t) [Data.due_dt_ex{1}(xq, yq, t), Data.due_dt_ex{2}(xq, yq, t)], ...
                                @(xq, yq, t, lam) -lam.*(Data.grad_ue_ex{1}(xq,yq, t) + Data.grad_ue_ex{4}(xq,yq,t))};
  
    SolutionInfo.t = t;
    SolutionInfo.label = 'E';

    [Xh, Xexact] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
 
end

