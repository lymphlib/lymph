%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 9 October 2024
%> @brief  Compute solution for post processing
%>
%==========================================================================
%> @section classPoroElaAcuGetSolutionPostProcessing Class description
%> @brief  Compute solution for post processing
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param Solutions  Struct with solution vectors
%> @param time       time instant
%
%> @retval Xh        Cell array with one cell for each physics
%>                   (1: porous medium; 2: acoustic medium; 3: elasticity) containing a
%>                   struct with approximate solutions at quad nodes and
%>                   sequence of strings useful for post-processing
%> @retval Xexact    Cell array with one cell for each physics
%>                   (1: porous medium; 2: acoustic medium; 3: elasticity) containing a
%>                   struct with exact solutions at quad nodes and
%>                   and sequence of strings useful for post-processing
%>
%==========================================================================

function [Xh, Xexact] = GetSolutionPostProcessing(Data, femregion, neighbor, Solutions, time)

%% Nodes values
ref_qNodes_2D = linspace(0,1,Data.NPtsVisualization+2)';
ref_qNodes_2D = ref_qNodes_2D(2:end-1);
ref_qNodes_2D = [repmat(ref_qNodes_2D,Data.NPtsVisualization,1), reshape(repmat(ref_qNodes_2D,1,Data.NPtsVisualization)',Data.NPtsVisualization^2,1)];

%% Setup
% Points coordinates
Xp  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
Yp  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
Xa  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
Ya  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
Xe  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
Ye  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

% Solutions of poroelastic problem
Uph  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
dUph = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
Wph  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
dWph = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
Pph  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

% Solutions of acoustic problem
Uah  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
dUah = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
Pah  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

% Solutions of elastic problem
Ueh  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
dUeh = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
Peh  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

if Data.PlotExact

    % Exact solutions of poroelastic problem
    Up_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    dUp_ex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    Wp_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    dWp_ex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    Pp_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

    % Exact solutions of acoustic problem
    Ua_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    dUa_ex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Pa_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

    % Exact solutions of elastic problem
    Ue_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    dUe_ex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    Pe_ex  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);

end

% Points properties
Idp = zeros(femregion.nel_p,2);
Bdp = cell(femregion.nel_p,1);

Ida = zeros(femregion.nel_a,2);
Bda = cell(femregion.nel_a,1);

Ide = zeros(femregion.nel_e,2);
Bde = cell(femregion.nel_e,1);

% shift index
kindp = 0;
kinda = 0;
kinde = 0;

id_shift = max(Data.TagElPoro);
if(isempty(id_shift)); id_shift = 0; end

id_shift_e = max([Data.TagElAcu,Data.TagElPoro]);
if(isempty(id_shift_e)); id_shift_e = 0; end


%% Loop over the elements
% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie = 1:femregion.nel % loop over elements
    
    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    % Element id
    id_ie  = femregion.id(ie);
    
    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    
    % Extraction of element geometrical information
    coords_ie = femregion.coords_element{ie};
 
    % Poroelastic element
    if femregion.tag(ie) == 'P'
        
        % Counter boundary nodes
        CountBD_min = 0;
        CountBD_max = 0;

        % First node of the current element
        Idp(ie,1) = kindp + 1;

        % Local solution
        u1_loc = Solutions.up_h(index);
        u2_loc = Solutions.up_h(index+femregion.ndof_p);
        w1_loc = Solutions.wp_h(index);
        w2_loc = Solutions.wp_h(index+femregion.ndof_p);
        % Local velocity
        du1_loc = Solutions.dot_up_h(index);
        du2_loc = Solutions.dot_up_h(index+femregion.ndof_p);
        dw1_loc = Solutions.dot_wp_h(index);
        dw2_loc = Solutions.dot_wp_h(index+femregion.ndof_p);
        
        % Computation of nodes on the physical element
        qNodes_2D = ref_qNodes_2D;

        qNodes_2D(:,1) = ref_qNodes_2D(:,1)*(femregion.bbox(ie,2)-femregion.bbox(ie,1))+femregion.bbox(ie,1);
        qNodes_2D(:,2) = ref_qNodes_2D(:,2)*(femregion.bbox(ie,4)-femregion.bbox(ie,3))+femregion.bbox(ie,3);
        
        if Data.NPtsVisualization > 1
            dist_2D = norm(qNodes_2D(end,:)-qNodes_2D(end-1,:));
        end

        qNodes_2D = qNodes_2D(inpolygon(qNodes_2D(:,1),qNodes_2D(:,2),coords_ie(:,1),coords_ie(:,2)),:);

        xq  = qNodes_2D(:,1);
        yq  = qNodes_2D(:,2);
        lqn = length(xq);

        % Construction and evaluation of the basis functions in the visualization points
        [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

        % Evaluation of physical parameters
        beta = Data.beta{id_ie}(xq,yq);
        m = Data.m{id_ie}(xq,yq);

        % Approximated solutions at quadrature points
        up1h_loc  = phiq*u1_loc;
        up2h_loc  = phiq*u2_loc;
        dup1h_loc  = phiq*du1_loc;
        dup2h_loc  = phiq*du2_loc;

        wp1h_loc  = phiq*w1_loc;
        wp2h_loc  = phiq*w2_loc;
        dwp1h_loc  = phiq*dw1_loc;
        dwp2h_loc  = phiq*dw2_loc;

        local_u_div = gradqx*u1_loc + gradqy*u2_loc;
        local_w_div = gradqx*w1_loc + gradqy*w2_loc;

        pp_loc =  - m .* (beta .* local_u_div + local_w_div);

        % Fill the output Ue, dUe and Pe
        Xp(kindp + 1 : kindp + lqn,:) = xq;
        Yp(kindp + 1 : kindp + lqn,:) = yq;
        Uph(kindp + 1 : kindp + lqn,:) = [up1h_loc, up2h_loc];
        Wph(kindp + 1 : kindp + lqn,:) = [wp1h_loc, wp2h_loc];
        Pph(kindp + 1 : kindp + lqn,:) = pp_loc;
        dUph(kindp + 1 : kindp + lqn,:) = [dup1h_loc, dup2h_loc];
        dWph(kindp + 1 : kindp + lqn,:) = [dwp1h_loc, dwp2h_loc];

        if Data.PlotExact
            % Exact solutions at quadrature points
            up1ex_loc = Data.up_ex{1}(xq,yq)*Data.up_t_ex{1}(time);
            up2ex_loc = Data.up_ex{2}(xq,yq)*Data.up_t_ex{1}(time);
            dup1ex_loc = Data.up_ex{1}(xq,yq)*Data.dup_t_ex{1}(time);
            dup2ex_loc = Data.up_ex{2}(xq,yq)*Data.dup_t_ex{1}(time);

            ppex_loc = - m .* (beta .* (Data.grad_up_ex{1}(xq,yq)+Data.grad_up_ex{2}(xq,yq))*Data.up_t_ex{1}(time) ...
                + (Data.grad_wp_ex{1}(xq,yq)+Data.grad_wp_ex{2}(xq,yq))*Data.wp_t_ex{1}(time));

            wp1ex_loc = Data.wp_ex{1}(xq,yq)*Data.wp_t_ex{1}(time);
            wp2ex_loc = Data.wp_ex{2}(xq,yq)*Data.wp_t_ex{1}(time);
            dwp1ex_loc = Data.wp_ex{1}(xq,yq)*Data.dwp_t_ex{1}(time);
            dwp2ex_loc = Data.wp_ex{2}(xq,yq)*Data.dwp_t_ex{1}(time);

            Up_ex(kindp + 1 : kindp + lqn,:) = [up1ex_loc, up2ex_loc];
            Wp_ex(kindp + 1 : kindp + lqn,:) = [wp1ex_loc, wp2ex_loc];
            Pp_ex(kindp + 1 : kindp + lqn,:) =  ppex_loc;
            dUp_ex(kindp + 1 : kindp + lqn,:) = [dup1ex_loc, dup2ex_loc];
            dWp_ex(kindp + 1 : kindp + lqn,:) = [dwp1ex_loc, dwp2ex_loc];

        end

        % Update auxiliary indexes
        kindp = kindp + lqn ;
        CountBD_min = CountBD_min + lqn + 1;
        CountBD_max = CountBD_max + lqn;

        % Loop over faces
        for iedg = 1 : neighbor.nedges(ie)

            % Extraction of the edge coordinates
            if iedg == neighbor.nedges(ie)
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(1,:);
            else
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(iedg+1,:);
            end

            if Data.NPtsVisualization > 1
                ref_qNodes_1D = linspace(0,1,ceil(norm(p1-p2)/dist_2D)+1)';
                ref_qNodes_1D = ref_qNodes_1D(1:end-1);

                [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);
            else
                qNodes_1D = p1;
            end

            % Construction of quadrature nodes on the face
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            lqn = length(xq);

            % Evaluation of physical parameters
            beta = Data.beta{id_ie}(xq,yq);
            m = Data.m{id_ie}(xq,yq);

            % Construction and evaluation of the basis functions in the visualization points
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

            % Approximated solutions at quadrature points
            up1h_loc  = phiedgeq*u1_loc;
            up2h_loc  = phiedgeq*u2_loc;
            dup1h_loc = phiedgeq*du1_loc;
            dup2h_loc = phiedgeq*du2_loc;

            wp1h_loc  = phiedgeq*w1_loc;
            wp2h_loc  = phiedgeq*w2_loc;
            dwp1h_loc = phiedgeq*dw1_loc;
            dwp2h_loc = phiedgeq*dw2_loc;

            local_u_div = gradedgeqx*u1_loc + gradedgeqy*u2_loc;
            local_w_div = gradedgeqx*w1_loc + gradedgeqy*w2_loc;

            pp_loc =  - m .* (beta .* local_u_div + local_w_div);

            % Fill the output Ue, dUe and Pe
            Xp(kindp + 1 : kindp + lqn,:) = xq;
            Yp(kindp + 1 : kindp + lqn,:) = yq;
            Uph(kindp + 1 : kindp + lqn,:) = [up1h_loc, up2h_loc];
            Wph(kindp + 1 : kindp + lqn,:) = [wp1h_loc, wp2h_loc];
            Pph(kindp + 1 : kindp + lqn,:) = pp_loc;
            dUph(kindp + 1 : kindp + lqn,:) = [dup1h_loc, dup2h_loc];
            dWph(kindp + 1 : kindp + lqn,:) = [dwp1h_loc, dwp2h_loc];

            if Data.PlotExact
                % Exact solutions at quadrature points
                up1ex_loc = Data.up_ex{1}(xq,yq)*Data.up_t_ex{1}(time);
                up2ex_loc = Data.up_ex{2}(xq,yq)*Data.up_t_ex{1}(time);
                dup1ex_loc = Data.up_ex{1}(xq,yq)*Data.dup_t_ex{1}(time);
                dup2ex_loc = Data.up_ex{2}(xq,yq)*Data.dup_t_ex{1}(time);

                ppex_loc = - m .* (beta .* (Data.grad_up_ex{1}(xq,yq)+Data.grad_up_ex{2}(xq,yq))*Data.up_t_ex{1}(time) ...
                    + (Data.grad_wp_ex{1}(xq,yq)+Data.grad_wp_ex{2}(xq,yq))*Data.wp_t_ex{1}(time));

                wp1ex_loc = Data.wp_ex{1}(xq,yq)*Data.wp_t_ex{1}(time);
                wp2ex_loc = Data.wp_ex{2}(xq,yq)*Data.wp_t_ex{1}(time);
                dwp1ex_loc = Data.wp_ex{1}(xq,yq)*Data.dwp_t_ex{1}(time);
                dwp2ex_loc = Data.wp_ex{2}(xq,yq)*Data.dwp_t_ex{1}(time);

                Up_ex(kindp + 1 : kindp + lqn,:) = [up1ex_loc, up2ex_loc];
                Wp_ex(kindp + 1 : kindp + lqn,:) = [wp1ex_loc, wp2ex_loc];
                Pp_ex(kindp + 1 : kindp + lqn,:) =  ppex_loc;
                dUp_ex(kindp + 1 : kindp + lqn,:) = [dup1ex_loc, dup2ex_loc];
                dWp_ex(kindp + 1 : kindp + lqn,:) = [dwp1ex_loc, dwp2ex_loc];

            end

            % Update auxiliary indexes
            kindp  = kindp  + lqn;
            CountBD_max = CountBD_max + lqn;
        end

        % save a vector from first to last (index of) boundary node of the current element
        Bdp{ie} = (CountBD_min:CountBD_max)';

        % Last node of the current element
        Idp(ie,2) = kindp;

    end
    
    % Elastic element
    if femregion.tag(ie) == 'A'
        
        % Counter boundary nodes
        CountBD_min = 0;
        CountBD_max = 0;

        ie_a = ie-femregion.nel_p;
        % First node of the current element
        Ida(ie_a,1) = kinda + 1;
   
        index_a = index - femregion.ndof_p;
        % Local solution
        phi_loc = Solutions.phi_h(index_a);
        % Local velocity
        dphi_loc = Solutions.dot_phi_h(index_a);
        
        % Computation of nodes on the physical element
        qNodes_2D = ref_qNodes_2D;

        qNodes_2D(:,1) = ref_qNodes_2D(:,1)*(femregion.bbox(ie,2)-femregion.bbox(ie,1))+femregion.bbox(ie,1);
        qNodes_2D(:,2) = ref_qNodes_2D(:,2)*(femregion.bbox(ie,4)-femregion.bbox(ie,3))+femregion.bbox(ie,3);
        
        if Data.NPtsVisualization > 1
            dist_2D = norm(qNodes_2D(end,:)-qNodes_2D(end-1,:));
        end

        qNodes_2D = qNodes_2D(inpolygon(qNodes_2D(:,1),qNodes_2D(:,2),coords_ie(:,1),coords_ie(:,2)),:);
   
        xq  = qNodes_2D(:,1);
        yq  = qNodes_2D(:,2);
        lqn = length(xq);
            
        % Construction and evaluation of the basis functions in the visualization points
        phiq = Evalshape2D(femregion, ie, qNodes_2D);

        % Evaluation of physical parameters
        rho_a = Data.rho_a{id_ie-id_shift}(xq,yq);

        % Approximated solutions at quadrature points
        phih_loc  = phiq*phi_loc;
        dphih_loc  = phiq*dphi_loc;

        Xa(kinda + 1 : kinda + lqn,:) = xq;
        Ya(kinda + 1 : kinda + lqn,:) = yq;

        % Fill the output Ua, dUa and Pa
        Uah(kinda + 1 : kinda + lqn,:) = phih_loc;
        Pah(kinda + 1 : kinda + lqn,:) = rho_a .* dphih_loc;
        dUah(kinda + 1 : kinda + lqn,:) = dphih_loc;

        if Data.PlotExact
            % Exact solutions at quadrature points
            phiex_loc  = Data.phi_ex{1}(xq,yq)*Data.phi_t_ex{1}(time);
            dphiex_loc = Data.phi_ex{1}(xq,yq)*Data.dphi_t_ex{1}(time);

            % Fill the output Ua, dUa and Pa
            Ua_ex(kinda + 1 : kinda + lqn,:) = phiex_loc;
            Pa_ex(kinda + 1 : kinda + lqn,:) = rho_a .* dphiex_loc;
            dUa_ex(kinda + 1 : kinda + lqn,:) = dphiex_loc;

        end

        kinda = kinda + lqn ;
        CountBD_min = CountBD_min + lqn + 1;
        CountBD_max = CountBD_max + lqn;

        % Loop over faces
        for iedg = 1 : neighbor.nedges(ie)

            % Extraction of the edge coordinates
            if iedg == neighbor.nedges(ie)
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(1,:);
            else
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(iedg+1,:);
            end

            if Data.NPtsVisualization > 1
                ref_qNodes_1D = linspace(0,1,ceil(norm(p1-p2)/dist_2D)+1)';
                ref_qNodes_1D = ref_qNodes_1D(1:end-1);

                [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);
            else
                qNodes_1D = p1;
            end

            % Construction of quadrature nodes on the face
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            lqn = length(xq);

            % Evaluation of physical parameters
            rho_a = Data.rho_a{id_ie-id_shift}(xq,yq);
            
            % Construction and evaluation of the basis functions in the visualization points
            phiedgeq = Evalshape2D(femregion, ie, qNodes_1D);

             % Approximated solutions at quadrature points
            phih_loc  = phiedgeq*phi_loc;
            dphih_loc  = phiedgeq*dphi_loc;
            
            Xa(kinda + 1 : kinda + lqn,:) = xq;
            Ya(kinda + 1 : kinda + lqn,:) = yq;
            
            % Fill the output Ua, dUa and Pa
            Uah(kinda + 1 : kinda + lqn,:) = phih_loc;
            Pah(kinda + 1 : kinda + lqn,:) = rho_a .* dphih_loc;
            dUah(kinda + 1 : kinda + lqn,:) = dphih_loc;
            
            if Data.PlotExact   
                % Exact solutions at quadrature points
                phiex_loc  = Data.phi_ex{1}(xq,yq)*Data.phi_t_ex{1}(time);
                dphiex_loc = Data.phi_ex{1}(xq,yq)*Data.dphi_t_ex{1}(time);

                % Fill the output Ua, dUa and Pa
                Ua_ex(kinda + 1 : kinda + lqn,:) = phiex_loc;
                Pa_ex(kinda + 1 : kinda + lqn,:) = rho_a .* dphiex_loc;
                dUa_ex(kinda + 1 : kinda + lqn,:) = dphiex_loc;

            end           

            kinda = kinda + lqn ;
            CountBD_max = CountBD_max + lqn;

        end

        % save a vector from first to last (index of) boundary node of the current element
        Bda{ie_a} = (CountBD_min:CountBD_max)';

        % Last node of the current element
        Ida(ie_a,2) = kinda;

    end
    
    
    % Elastic element
    if femregion.tag(ie) == 'E'
        
        % Counter boundary nodes
        CountBD_min = 0;
        CountBD_max = 0;

        ie_e = ie - femregion.nel_p - femregion.nel_a;

        % First node of the current element
        Ide(ie_e,1) = kinde + 1;
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        % Local solution
        u1_loc = Solutions.ue_h(index_e);
        u2_loc = Solutions.ue_h(index_e+femregion.ndof_e);
        % Local velocity
        du1_loc = Solutions.dot_ue_h(index_e);
        du2_loc = Solutions.dot_ue_h(index_e+femregion.ndof_e);
        
        % Computation of nodes on the physical element
        qNodes_2D = ref_qNodes_2D;

        qNodes_2D(:,1) = ref_qNodes_2D(:,1)*(femregion.bbox(ie,2)-femregion.bbox(ie,1))+femregion.bbox(ie,1);
        qNodes_2D(:,2) = ref_qNodes_2D(:,2)*(femregion.bbox(ie,4)-femregion.bbox(ie,3))+femregion.bbox(ie,3);
        
        if Data.NPtsVisualization > 1
            dist_2D = norm(qNodes_2D(end,:)-qNodes_2D(end-1,:));
        end

        qNodes_2D = qNodes_2D(inpolygon(qNodes_2D(:,1),qNodes_2D(:,2),coords_ie(:,1),coords_ie(:,2)),:);
    
        xq  = qNodes_2D(:,1);
        yq  = qNodes_2D(:,2);
        lqn = length(xq);

        % Construction and evaluation of the basis functions in the visualization points
        [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

        % Evaluation of physical parameters
        lambda = Data.lam_el{id_ie-id_shift_e}(xq,yq);

        % Approximated solutions at quadrature points
        ue1h_loc  = phiq*u1_loc;
        ue2h_loc  = phiq*u2_loc;
        pe_loc    = - lambda .* (gradqx*u1_loc + gradqy*u2_loc);

        due1h_loc  = phiq*du1_loc;
        due2h_loc  = phiq*du2_loc;

        Xe(kinde + 1 : kinde + lqn,:) = xq;
        Ye(kinde + 1 : kinde + lqn,:) = yq;

        % Fill the output Ue, dUe and Pe
        Ueh(kinde + 1 : kinde + lqn,:)  = [ue1h_loc, ue2h_loc];
        Peh(kinde + 1 : kinde + lqn,:)  = pe_loc;
        dUeh(kinde + 1 : kinde + lqn,:) = [due1h_loc, due2h_loc];

        if Data.PlotExact
            % Exact solutions at quadrature points
            ue1ex_loc = Data.ue_ex{1}(xq,yq)*Data.ue_t_ex{1}(time);
            ue2ex_loc = Data.ue_ex{2}(xq,yq)*Data.ue_t_ex{1}(time);
            due1ex_loc = Data.ue_ex{1}(xq,yq)*Data.due_t_ex{1}(time);
            due2ex_loc = Data.ue_ex{2}(xq,yq)*Data.due_t_ex{1}(time);

            peex_loc = - lambda .* (Data.grad_ue_ex{1}(xq,yq)+Data.grad_ue_ex{2}(xq,yq))*Data.ue_t_ex{1}(time);

            Ue_ex(kinde + 1 : kinde + lqn,:)  = [ue1ex_loc, ue2ex_loc];
            Pe_ex(kinde + 1 : kinde + lqn,:)  = peex_loc;
            dUe_ex(kinde + 1 : kinde + lqn,:) = [due1ex_loc, due2ex_loc];

        end

        kinde = kinde + lqn;
        CountBD_min = CountBD_min + lqn + 1;
        CountBD_max = CountBD_max + lqn;

        % Loop over faces
        for iedg = 1 : neighbor.nedges(ie)

            % Extraction of the edge coordinates
            if iedg == neighbor.nedges(ie)
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(1,:);
            else
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(iedg+1,:);
            end

            if Data.NPtsVisualization > 1
                ref_qNodes_1D = linspace(0,1,ceil(norm(p1-p2)/dist_2D)+1)';
                ref_qNodes_1D = ref_qNodes_1D(1:end-1);

                [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);
            else
                qNodes_1D = p1;
            end

            % Construction of quadrature nodes on the face
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            lqn = length(xq);

            % Construction and evaluation of the basis functions in the visualization points
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

            % Evaluation of physical parameters
            lambda = Data.lam_el{id_ie-id_shift_e}(xq,yq);
           
            % Approximated solutions at quadrature points
            ue1h_loc  = phiedgeq*u1_loc;
            ue2h_loc  = phiedgeq*u2_loc;
            pe_loc    = - lambda .* (gradedgeqx*u1_loc + gradedgeqy*u2_loc);
            
            due1h_loc  = phiedgeq*du1_loc;
            due2h_loc  = phiedgeq*du2_loc;
           
            Xe(kinde + 1 : kinde + lqn,:) = xq;
            Ye(kinde + 1 : kinde + lqn,:) = yq;
            
            % Fill the output Ue, dUe and Pe
            Ueh(kinde + 1 : kinde + lqn,:)  = [ue1h_loc, ue2h_loc];
            Peh(kinde + 1 : kinde + lqn,:)  = pe_loc;
            dUeh(kinde + 1 : kinde + lqn,:) = [due1h_loc, due2h_loc];
            
            if Data.PlotExact
                % Exact solutions at quadrature points
                ue1ex_loc = Data.ue_ex{1}(xq,yq)*Data.ue_t_ex{1}(time);
                ue2ex_loc = Data.ue_ex{2}(xq,yq)*Data.ue_t_ex{1}(time);
                due1ex_loc = Data.ue_ex{1}(xq,yq)*Data.due_t_ex{1}(time);
                due2ex_loc = Data.ue_ex{2}(xq,yq)*Data.due_t_ex{1}(time);

                peex_loc = - lambda .* (Data.grad_ue_ex{1}(xq,yq)+Data.grad_ue_ex{2}(xq,yq))*Data.ue_t_ex{1}(time);

                Ue_ex(kinde + 1 : kinde + lqn,:)  = [ue1ex_loc, ue2ex_loc];
                Pe_ex(kinde + 1 : kinde + lqn,:)  = peex_loc;
                dUe_ex(kinde + 1 : kinde + lqn,:) = [due1ex_loc, due2ex_loc];
            
            end        

            kinde = kinde + lqn ;
            CountBD_max = CountBD_max + lqn;

        end

        % save a vector from first to last (index of) boundary node of the current element
        Bde{ie_e} = (CountBD_min:CountBD_max)';

        % Last node of the current element
        Ide(ie_e,2) = kinde;
    end
    
    
    
end

%% Output
Xh = {};

XPlotp{1} = Xp(1:kindp);
XPlotp{2} = Yp(1:kindp);
XPlotp{3} = Uph(1:kindp,:);
XPlotp{4} = Wph(1:kindp,:);
XPlotp{5} = Pph(1:kindp,:);
XPlotp{6} = dUph(1:kindp,:);
XPlotp{7} = dWph(1:kindp,:);
Xh{1}.Solution = XPlotp;
Xh{1}.StrVTK  = {'x','y','uph','wph','pph','duph','dwph'};
Xh{1}.StrPlot = {'x','y','up','wp','pp','dup','dwp'};
Xh{1}.StrCSV = {'x','y','uph1','uph2','wph1','wph2','pph','duph1','duph2','dwph1','dwph2'};
Xh{1}.Id       = Idp;
Xh{1}.Bd       = Bdp;

XPlota{1} = Xa(1:kinda);
XPlota{2} = Ya(1:kinda);
XPlota{3} = Uah(1:kinda,:);
XPlota{4} = Pah(1:kinda,:);
XPlota{5} = dUah(1:kinda,:);
Xh{2}.Solution = XPlota;
Xh{2}.StrVTK  = {'x','y','phiah','pah','dphiah'};
Xh{2}.StrPlot = {'x','y','\varphi_a','pa','\varphi_{a,t}'};
Xh{2}.StrCSV = {'x','y','phiah','pah','dphiah'};
Xh{2}.Id       = Ida;
Xh{2}.Bd       = Bda;

XPlote{1} = Xe(1:kinde);
XPlote{2} = Ye(1:kinde);
XPlote{3} = Ueh(1:kinde,:);
XPlote{4} = Peh(1:kinde,:);
XPlote{5} = dUeh(1:kinde,:);
Xh{3}.Solution = XPlote;
Xh{3}.StrVTK  = {'x','y','ueh','peh','veh'};
Xh{3}.StrPlot = {'x','y','ue','pe','ve'};
Xh{3}.StrCSV = {'x','y','ueh1','ueh2','peh','veh1','veh2'};
Xh{3}.Id       = Ide;
Xh{3}.Bd       = Bde;

if Data.PlotExact
    Xexact = {};

    XPlotExactp{1} = Xp(1:kindp);
    XPlotExactp{2} = Yp(1:kindp);
    XPlotExactp{3} = Up_ex(1:kindp,:);
    XPlotExactp{4} = Wp_ex(1:kindp,:);
    XPlotExactp{5} = Pp_ex(1:kindp,:);
    XPlotExactp{6} = dUp_ex(1:kindp,:);
    XPlotExactp{7} = dWp_ex(1:kindp,:);
    Xexact{1}.Solution = XPlotExactp;
    Xexact{1}.StrPlot  = Xh{1}.StrPlot;

    XPlotExacta{1} = Xa(1:kinda);
    XPlotExacta{2} = Ya(1:kinda);
    XPlotExacta{3} = Ua_ex(1:kinda,:);
    XPlotExacta{4} = Pa_ex(1:kinda,:);
    XPlotExacta{5} = dUa_ex(1:kinda,:);
    Xexact{2}.Solution = XPlotExacta;
    Xexact{2}.StrPlot  = Xh{2}.StrPlot;

    XPlotExacte{1} = Xe(1:kinde);
    XPlotExacte{2} = Ye(1:kinde);
    XPlotExacte{3} = Ue_ex(1:kinde,:);
    XPlotExacte{4} = Pe_ex(1:kinde,:);
    XPlotExacte{5} = dUe_ex(1:kinde,:);
    Xexact{3}.Solution = XPlotExacte;
    Xexact{3}.StrPlot  = Xh{3}.StrPlot;

else
    Xexact = {[],[],[]};
end


