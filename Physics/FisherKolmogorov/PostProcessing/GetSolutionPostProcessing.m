%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti, Caterina Leimer Saglio
%> @date 5 February 2026
%> @brief Compute the solution for post-processing for the heat equation.
%>
%==========================================================================
%> @section classFKPPGetSolutionPostProcessing Class description
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

    SolutionInfo.outputnames = {'c','p'};
    SolutionInfo.outputsizes = { 1 , 1 };
    SolutionInfo.outputsol   = { Solution.c_h , femregion.degree };
    SolutionInfo.outputexpr  = { @(phi, gradx, grady, c) phi*c        , @(phi, gradx, grady, p) p };
    SolutionInfo.outputexact = { @(xq, yq, t) Data.c_ex{1}(xq, yq, t) , @(xq, yq, t) 0*xq };

    if Data.Adaptivity
        SolutionInfo.outputnames(3:7) = {'tau','tauE','tauJ','tauN','tauT'};
        SolutionInfo.outputsizes(3:7) = {  1  ,  1   ,  1   ,  1   ,  1   };
        SolutionInfo.outputsol(3:7)   = { Solution.Indicator.tau , ...
                                          Solution.Indicator.tau_E, ...
                                          Solution.Indicator.tau_J, ...
                                          Solution.Indicator.tau_N, ...
                                          Solution.Indicator.tau_T};

        SolutionInfo.outputexpr(3:7)  = {@(phi, gradx, grady, tau) tau};
        SolutionInfo.outputexact(3:7) = {@(xq, yq, t) 0*xq};
    end
  
    SolutionInfo.t = t;
    SolutionInfo.label = Data.LabEl{1};


    [Xh, Xexact] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
 
end
