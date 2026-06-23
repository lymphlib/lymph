%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti, Caterina Leimer Saglio
%> @date 28 October 2025
%> @brief Compute the solution for post-processing for the Laplace equation.
%>
%==========================================================================
%> @section classLaplacianGetSolutionPostProcessing Class description
%==========================================================================
%> @brief  Compute the solution for post-processing for the Laplace equation.
%>
%> @param Data        Struct with problem's data
%> @param femregion   Struct for finite elements data
%> @param neighbor    Struct for neighboring elements data
%> @param Solution    Modal solution and adaptive indicators
%>
%> @retval Xh         Cell array containing a struct with approximate
%>                    solutions strings useful for post processing
%> @retval Xexact     Cell array containing a struct with exact
%>                    solutions strings useful for post processing
%>
%==========================================================================

function [Xh, Xexact] = GetSolutionPostProcessing(Data, femregion, neighbor, Solution)

    SolutionInfo.outputnames = {'u','p'};
    SolutionInfo.outputsizes = { 1 , 1 };
    SolutionInfo.outputsol   = { Solution.U , femregion.degree };
    SolutionInfo.outputexpr  = { @(phi, gradx, grady, u) phi*u    , @(phi, gradx, grady, p) p };
    SolutionInfo.outputexact = { @(xq, yq, t) Data.u_ex{1}(xq,yq) , @(xq, yq, t) 0*xq };

    if Data.Adaptivity
        newNames  = {'tau','tauE','tauJ','tauN','tauT'};
        newSizes  = { 1,    1,     1,     1,     1 };
        newSols   = { ...
            Solution.Indicator.tau, ...
            Solution.Indicator.tau_E, ...
            Solution.Indicator.tau_J, ...
            Solution.Indicator.tau_N, ...
            Solution.Indicator.tau_T ...
        };
    
        SolutionInfo.outputnames = [SolutionInfo.outputnames, newNames];
        SolutionInfo.outputsizes = [SolutionInfo.outputsizes, newSizes];
        SolutionInfo.outputsol   = [SolutionInfo.outputsol,   newSols];
    
        for k = 1:numel(newNames)
            SolutionInfo.outputexpr{end+1}  = @(phi,gradx,grady,val) val;
            SolutionInfo.outputexact{end+1} = @(xq,yq,t) 0*xq;
        end
    end
  
    % Fictitious time
    SolutionInfo.t = 0;
    SolutionInfo.label = Data.LabEl{1};

    [Xh, Xexact] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
 
end
