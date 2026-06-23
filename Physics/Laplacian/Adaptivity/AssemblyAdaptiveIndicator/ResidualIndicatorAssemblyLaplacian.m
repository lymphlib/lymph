%> @file   ResidualIndicatorAssemblyLaplacian.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   12 May 2026
%> @brief Computation of the residual indicator for Laplacian problem.
%>
%==========================================================================
%> @section classResidualIndicatorAssemblyLaplacian description
%==========================================================================
%> @brief            Computation of the residual indicator for Laplacian problem.
%>
%> @param Data       Struct with problem's data
%> @param Indicator  Indicator struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Indicator Indicator struct
%>                   
%==========================================================================

function [Indicator] = ResidualIndicatorAssemblyLaplacian(Data, Indicator, elem, ie, id, nbases, AssembInfo)
        
        % Physical parameters evaluation    
        mu      = Data.mu{1}(elem.xq,elem.yq);
        f_loc   = Data.source{1}(elem.xq,elem.yq);

        lap_loc = (mu .* (elem.lapqxx + elem.lapqyy)) * AssembInfo.uh(elem.index);

        % Local residual integral assembly
        Indicator.tau_E = Data.h.^2 * (elem.dx'* (f_loc + lap_loc).^2);
        
end
                    
