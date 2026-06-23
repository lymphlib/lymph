%> @file   ResidualIndicatorAssemblyFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief Computation of the residual indicator for FKPP equation.
%>
%==========================================================================
%> @section classResidualIndicatorAssemblyFKPP description
%==========================================================================
%> @brief            Computation of the residual indicator for FKPP equation.
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

function [Indicator] = ResidualIndicatorAssemblyFKPP(Data, Indicator, elem, ie, id, nbases, AssembInfo)
        
        % Physical parameters evaluation    
        Dext    = Data.D_ext{1}(elem.xq,elem.yq);
        alpha   = Data.alpha{1}(elem.xq,elem.yq);
        f_loc   = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t,Dext,alpha);
        f_old_loc   = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t-Data.dt,Dext,alpha);

        lap_loc     = (Dext .* (elem.lapqxx + elem.lapqyy)) * AssembInfo.c_h(elem.index);
        lap_old_loc = (Dext .* (elem.lapqxx + elem.lapqyy)) * AssembInfo.c_old(elem.index);
        
        dtu_loc = (elem.phiq*AssembInfo.c_h(elem.index)-elem.phiq*AssembInfo.c_old(elem.index))./Data.dt;

        reac_loc      = ((alpha .* elem.phiq) * AssembInfo.c_h(elem.index)).*(1 - elem.phiq * AssembInfo.c_h(elem.index));
        reac_old_loc  = ((alpha .* elem.phiq) * AssembInfo.c_old(elem.index)).*(1 - elem.phiq * AssembInfo.c_old(elem.index));
        reac_mix1_loc = ((alpha .* elem.phiq) * AssembInfo.c_h(elem.index)).*(1 - elem.phiq * AssembInfo.c_old(elem.index));
        reac_mix2_loc = ((alpha .* elem.phiq) * AssembInfo.c_old(elem.index)).*(1 - elem.phiq * AssembInfo.c_h(elem.index));

        % Local residual integral assembly
        Indicator.tau_E = Data.h.^2 * (elem.dx'* (Data.theta*f_loc + (1-Data.theta)*f_old_loc ...
                                                    + Data.theta*lap_loc + (1-Data.theta)*lap_old_loc ...
                                                    + Data.theta^2*reac_loc + (1-Data.theta)^2*reac_old_loc ...
                                                    + Data.theta*(1-Data.theta)*(reac_mix1_loc+reac_mix2_loc) ...
                                                    - dtu_loc).^2);

end
                    
