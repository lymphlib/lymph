%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Compute the solution for post-processing for the heat equation.
%>
%==========================================================================
%> @section classFHNGetSolutionPostProcessing Class description
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

    SolutionInfo.outputnames = {'u','w','p'};
    SolutionInfo.outputsizes = { 1 ,1, 1 };
    SolutionInfo.outputsol   = { Solution.u_h, Solution.w_h, femregion.degree };
    SolutionInfo.outputexpr  = { @(phi, gradx, grady, c) phi*c, @(phi, gradx, grady, c) phi*c , @(phi, gradx, grady, p) p };
    SolutionInfo.outputexact = { @(xq, yq, t) Data.u_ex{1}(xq, yq, t),  @(xq, yq, t) Data.w_ex{1}(xq, yq, t) , @(xq, yq, t) 0*xq };

    if Data.Adaptivity
        SolutionInfo.outputnames(4:8) = {'tau','tauE','tauJ','tauN','tauT'};
        SolutionInfo.outputsizes(4:8) = {  1  ,  1   ,  1   ,  1   ,  1   };
        SolutionInfo.outputsol(4:8)   = { Solution.Indicator.tau , ...
                                          Solution.Indicator.tau_E, ...
                                          Solution.Indicator.tau_J, ...
                                          Solution.Indicator.tau_N, ...
                                          Solution.Indicator.tau_T};

        SolutionInfo.outputexpr(4:8)  = {@(phi, gradx, grady, tau) tau};
        SolutionInfo.outputexact(4:8) = {@(xq, yq, t) 0*xq};
    end
  
    SolutionInfo.t = t;
    SolutionInfo.label = Data.LabEl{1};
    
    [Xh, Xexact] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
 
end