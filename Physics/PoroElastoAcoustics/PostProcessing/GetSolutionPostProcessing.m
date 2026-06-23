%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti, Caterina Leimer Saglio
%> @date 5 June 2026
%> @brief Compute the solution for post-processing for the elastodynamics equation.
%>
%==========================================================================
%> @section classPoroAcuElaGetSolutionPostProcessing Class description
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


    %% Visualization of Poroelasticity
    if any(femregion.phys == 'P')
        id_phys = (femregion.phys == 'P');
        SolutionInfo.outputnames = {'u_p','du_p','w_p','dw_p','p_p'};
        SolutionInfo.outputsizes = { 2 , 2 , 2 , 2 , 1 };
        SolutionInfo.phys_param  = {@(xq, yq, t, id) Data.m{id}(xq,yq).*Data.beta{id}(xq,yq)};
        SolutionInfo.outputsol   = { [Solution.up_h(1:femregion.ndof_phys(id_phys)),     Solution.up_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))], ...
                                     [Solution.dot_up_h(1:femregion.ndof_phys(id_phys)), Solution.dot_up_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))], ...
                                     [Solution.wp_h(1:femregion.ndof_phys(id_phys)),     Solution.wp_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))], ...
                                     [Solution.dot_wp_h(1:femregion.ndof_phys(id_phys)), Solution.dot_wp_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))], ...
                                     [Solution.up_h(1:femregion.ndof_phys(id_phys)),     Solution.up_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys)), ...
                                      Solution.wp_h(1:femregion.ndof_phys(id_phys)),     Solution.wp_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))]};
        SolutionInfo.outputexpr  = {@(phi, gradx, grady, u)  [phi*u(:,1),  phi*u(:,2)], ...
                                    @(phi, gradx, grady, du) [phi*du(:,1), phi*du(:,2)], ...
                                    @(phi, gradx, grady, w)  [phi*w(:,1),  phi*w(:,2)], ...
                                    @(phi, gradx, grady, dw) [phi*dw(:,1), phi*dw(:,2)], ...
                                    @(phi, gradx, grady, u, par) -par.*(gradx*u(:,1)+grady*u(:,2)+gradx*u(:,3)+grady*u(:,4))};
        SolutionInfo.outputexact = {@(xq, yq, t) [Data.up_ex{1}(xq,yq)*Data.up_t_ex{1}(t),  Data.up_ex{2}(xq,yq)*Data.up_t_ex{1}(t)], ...
                                    @(xq, yq, t) [Data.up_ex{1}(xq,yq)*Data.dup_t_ex{1}(t), Data.up_ex{2}(xq,yq)*Data.dup_t_ex{1}(t)], ...
                                    @(xq, yq, t) [Data.wp_ex{1}(xq,yq)*Data.wp_t_ex{1}(t),  Data.wp_ex{2}(xq,yq)*Data.wp_t_ex{1}(t)], ...
                                    @(xq, yq, t) [Data.wp_ex{1}(xq,yq)*Data.dwp_t_ex{1}(t), Data.wp_ex{2}(xq,yq)*Data.dwp_t_ex{1}(t)], ...
                                    @(xq, yq, t, par) -par.*((Data.grad_up_ex{1}(xq,yq)+Data.grad_up_ex{2}(xq,yq))*Data.up_t_ex{1}(t)+(Data.grad_wp_ex{1}(xq,yq)+Data.grad_wp_ex{2}(xq,yq))*Data.wp_t_ex{1}(t))};
      
        SolutionInfo.t = t;
        SolutionInfo.label = 'P';
    
        [XhPoro, XexactPoro] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
    else
        XhPoro{1}     = {};
        XexactPoro{1} = {};
    end

    %% Visualization of Elastodynamics
    if any(femregion.phys == 'E')
        id_phys = (femregion.phys == 'E');
        SolutionInfo.outputnames = {'u_e','v_e','p_e'};
        SolutionInfo.outputsizes = { 2 , 2 , 1 };
        SolutionInfo.phys_param  = {@(xq, yq, t, id) Data.lam_el{id}(xq,yq)};
        SolutionInfo.outputsol   = { [Solution.ue_h(1:femregion.ndof_phys(id_phys)),     Solution.ue_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))], ...
                                     [Solution.dot_ue_h(1:femregion.ndof_phys(id_phys)), Solution.dot_ue_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))], ...
                                     [Solution.ue_h(1:femregion.ndof_phys(id_phys)),     Solution.ue_h(femregion.ndof_phys(id_phys)+1:2*femregion.ndof_phys(id_phys))]};
        SolutionInfo.outputexpr  = {@(phi, gradx, grady, u) [phi*u(:,1), phi*u(:,2)], ...
                                    @(phi, gradx, grady, v) [phi*v(:,1), phi*v(:,2)], ...
                                    @(phi, gradx, grady, u, lam) -lam.*(gradx*u(:,1)+grady*u(:,2))};
        SolutionInfo.outputexact = {@(xq, yq, t) [Data.ue_ex{1}(xq, yq)*Data.ue_t_ex{1}(t),  Data.ue_ex{2}(xq, yq)*Data.ue_t_ex{1}(t)], ...
                                    @(xq, yq, t) [Data.ue_ex{1}(xq, yq)*Data.due_t_ex{1}(t), Data.ue_ex{2}(xq, yq)*Data.due_t_ex{1}(t)], ...
                                    @(xq, yq, t, lam) -lam.*(Data.grad_ue_ex{1}(xq, yq) + Data.grad_ue_ex{4}(xq,yq))*Data.ue_t_ex{1}(t)};
      
        SolutionInfo.t = t;
        SolutionInfo.label = 'E';
    
        [XhEla, XexactEla] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
    else
        XhEla{1}     = {};
        XexactEla{1} = {};
    end

    %% Visualization of Acoustics
    if any(femregion.phys == 'A')
        id_phys = (femregion.phys == 'A');
        SolutionInfo.outputnames = {'\varphi_a', 'p_a', 'd\varphi_a'};
        SolutionInfo.outputsizes = { 1 , 1 , 1 };
        SolutionInfo.phys_param  = {@(xq, yq, t, id) Data.m{id}(xq,yq).*Data.beta{id}(xq,yq)};
        SolutionInfo.outputsol   = { Solution.phi_h(1:femregion.ndof_phys(id_phys)), ...
                                     Solution.dot_phi_h(1:femregion.ndof_phys(id_phys)), ...
                                     Solution.dot_phi_h(1:femregion.ndof_phys(id_phys))};
        SolutionInfo.outputexpr  = {@(phi, gradx, grady, p)  phi*p(:,1), ...
                                    @(phi, gradx, grady, dp, rho) rho.*(phi*dp(:,1)), ...
                                    @(phi, gradx, grady, dp) phi*dp(:,1)};
        SolutionInfo.outputexact = {@(xq, yq, t) Data.phi_ex{1}(xq,yq)*Data.phi_t_ex{1}(t), ...
                                    @(xq, yq, t, rho) rho.*Data.phi_ex{1}(xq,yq)*Data.dphi_t_ex{1}(t), ...
                                    @(xq, yq, t) Data.phi_ex{1}(xq,yq)*Data.dphi_t_ex{1}(t)};

        SolutionInfo.t = t;
        SolutionInfo.label = 'A';

        [XhAcu, XexactAcu] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo);
    else
        XhAcu{1}     = {};
        XexactAcu{1} = {};
    end

    Xh = {XhEla{1}, XhPoro{1}, XhAcu{1}};
    Xh = Xh(~cellfun('isempty', Xh));
    if Data.PlotExact
        Xexact = {XexactEla{1}, XexactPoro{1}, XexactAcu{1}};
        Xexact = Xexact(~cellfun('isempty', Xexact));
    else
        Xexact = {{}, {}, {}};
    end
end

