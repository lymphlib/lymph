%> @file  ComputeModalSolutionPoroAcuEla.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Compute modal expansions for poro-acoustic-elastic solution
%>
%==========================================================================
%> @section classComputeModalSolutionPoroAcuEla Class description
%> @brief  Compute Modal Solution for the Wave problem
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%
%> @retval [up,wp,phi_a,ue]  Vector with modal expansions
%>
%==========================================================================

function [up,wp,phi_a,ue] = ComputeModalSolutionPoroAcuEla(Data,femregion)
%% Quadrature values
[~, ~, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization
up1 = spalloc(femregion.ndof_p,1,femregion.ndof_p);
up2 = spalloc(femregion.ndof_p,1,femregion.ndof_p);
wp1 = spalloc(femregion.ndof_p,1,femregion.ndof_p);
wp2 = spalloc(femregion.ndof_p,1,femregion.ndof_p);

phi_a = spalloc(femregion.ndof_a,1,femregion.ndof_a);

ue1 = spalloc(femregion.ndof_e,1,femregion.ndof_e);
ue2 = spalloc(femregion.ndof_e,1,femregion.ndof_e);

%% Loop over the elements
% Visualization of computational progress
index_shift=0;
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie=1:femregion.nel
    
    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    % Element tag
    tag_ie = femregion.tag(ie);
    
    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    index_element = index_shift + (1:1:femregion.nedges(ie))';
    index_shift   = index_element(end);
    
    % Extraction of element geometrical information
    coords_ie          = femregion.coords_element{ie};
    
    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
    
    % Check if the element is poroelastic
    if tag_ie == 'P'

        % Inizialize local solution vector
        up1_loc = 0;
        up2_loc = 0;
        wp1_loc = 0;
        wp2_loc = 0;
    
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Construction of the basis functions
            [phiq, ~, ~] = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Compute exact solution
            up1ex_loc = Data.up_ex{1}(xq,yq);
            up2ex_loc = Data.up_ex{2}(xq,yq);
            wp1ex_loc = Data.wp_ex{1}(xq,yq);
            wp2ex_loc = Data.wp_ex{2}(xq,yq);
            
            up1_loc = up1_loc + (dx.*phiq)'*up1ex_loc;
            up2_loc = up2_loc + (dx.*phiq)'*up2ex_loc;
            wp1_loc = wp1_loc + (dx.*phiq)'*wp1ex_loc;
            wp2_loc = wp2_loc + (dx.*phiq)'*wp2ex_loc;
            
        end
        up1(index) = up1_loc;
        up2(index) = up2_loc;
        wp1(index) = wp1_loc;
        wp2(index) = wp2_loc;

    end

    % Check if the element is acoustic
    if tag_ie == 'A'
        
        index_a = index - femregion.ndof_p;
        % Inizialize local solution vector
        phi_a_loc = 0;
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Construction of the basis functions
            [phiq, ~, ~] = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Compute exact solution
            phiex_loc = Data.phi_ex{1}(xq,yq);
            
            phi_a_loc = phi_a_loc + (dx.*phiq)'*phiex_loc;
            
        end
        
        phi_a(index_a) = phi_a_loc;
        
    end

    % Check if the element is elastic
    if tag_ie == 'E'
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        % Inizialize local solution vector
        ue1_loc = 0;
        ue2_loc = 0;
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Construction of the basis functions
            [phiq, ~, ~] = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Compute exact solution
            ue1ex_loc = Data.ue_ex{1}(xq,yq);
            ue2ex_loc = Data.ue_ex{2}(xq,yq);
            
            ue1_loc = ue1_loc + (dx.*phiq)'*ue1ex_loc;
            ue2_loc = ue2_loc + (dx.*phiq)'*ue2ex_loc;
            
        end
        
        ue1(index_e) = ue1_loc;
        ue2(index_e) = ue2_loc;
        
    end
end


up = [up1; up2];
wp = [wp1; wp2];
ue = [ue1; ue2];

fprintf('\n');

