%> @file   ResidualIndicatorAssemblyFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief Computation of the residual indicator for FHN equation.
%>
%==========================================================================
%> @section classResidualIndicatorAssemblyFHN description
%==========================================================================
%> @brief            Computation of the residual indicator for FHN equation.
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
%>                    - \f$\tau_E\f$ = \f$\|h \theta(f + \nabla \cdot (\mu \nabla c_h) + (1 -\theta)(f^{old} + \nabla \cdot (\mu \nabla c_h^{old}) \|^2_{L2} \f$
%>                   
%==========================================================================

function [Indicator] = ResidualIndicatorAssemblyFHN(Data, Indicator, elem, ie, id, nbases, AssembInfo)
        
        % Physical parameters evaluation    
        Dext    = Data.D_ext{1}(elem.xq,elem.yq);
        a       = Data.a{1}(elem.xq,elem.yq);
        
        f_loc   = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t);
        f_old_loc   = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t-Data.dt);

        lap_loc     = (Dext .* (elem.lapqxx + elem.lapqyy)) * AssembInfo.u_h(elem.index);
        lap_old_loc = (Dext .* (elem.lapqxx + elem.lapqyy)) * AssembInfo.u_old(elem.index);
        
        dtu_loc = (elem.phiq*AssembInfo.u_h(elem.index)-elem.phiq*AssembInfo.u_old(elem.index))./Data.dt;

        reac_loc      = ((elem.phiq) * AssembInfo.u_h(elem.index)).*(elem.phiq * AssembInfo.u_h(elem.index) - 1).*(elem.phiq * AssembInfo.u_h(elem.index) - a) +  (elem.phiq * AssembInfo.w_h(elem.index) );
        reac_old_loc  = ((elem.phiq) * AssembInfo.u_old(elem.index)).*(elem.phiq * AssembInfo.u_old(elem.index) - 1).*(elem.phiq * AssembInfo.u_old(elem.index) - a) +  (elem.phiq * AssembInfo.w_old(elem.index) );

        % Local residual integral assembly
        Indicator.tau_E = Data.h.^2 * (elem.dx'* (Data.theta*f_loc + (1-Data.theta)*f_old_loc ...
                                                    + Data.theta*lap_loc + (1-Data.theta)*lap_old_loc ...
                                                    - Data.theta*reac_loc - (1-Data.theta)*reac_old_loc ...
                                                    - dtu_loc).^2);

end
                    
