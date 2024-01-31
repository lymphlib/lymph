%> @file  GetSolutionQuadPoints.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Compute solution at quadrature nodes
%>
%==========================================================================
%> @section classElastodynamicsGetSolutionQuadPoints Class description
%> @brief  Compute solution at quadrature nodes
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Solutions  Struct with solution vectors
%> @param time       time instant
%
%> @retval Displacement   Vector with solution displacement at quad nodes
%> @retval Pressure       Vector with solution pressure at quad nodes
%> @retval Velocity       Vector with solution velocity at quad nodes
%>
%==========================================================================

function [Displacement, Pressure, Velocity] = GetSolutionQuadPoints(Data, femregion, Solutions, time)

%% Quadrature values
[~, ~, ref_qNodes_2D, ~] = Quadrature(femregion.nqn);

%% Setup
% G = [ x | y | uh ... | uex ... ];
% Solid elastic displacement
Ue   = zeros(femregion.nel_e*femregion.nqn.^2,6);
% Solid elastic velcoty
dUe  = zeros(femregion.nel_e*femregion.nqn.^2,6);
% Pressure elastic
Pe   = zeros(femregion.nel_e*femregion.nqn.^2,3);

% shift index
kinde = 0;

id_shift_e = 0;


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
    
    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
    
    % Poroelastic element
    % Elastic element
    if femregion.tag(ie) == 'E'
        
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        % Local solution
        u1_loc = Solutions.ue_h(index_e);
        u2_loc = Solutions.ue_h(index_e+femregion.ndof_e);
        % Local velocity
        du1_loc = Solutions.dot_ue_h(index_e);
        du2_loc = Solutions.dot_ue_h(index_e+femregion.ndof_e);
        
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [~, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            lqn = length(xq);
            
            % Construction of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Evaluation of physical parameters
            lambda = Data.lam_el{id_ie-id_shift_e}(xq,yq);
            
            
            % Approximated solutions at quadrature points
            ue1h_loc  = phiq*u1_loc;
            ue2h_loc  = phiq*u2_loc;
            pe_loc    = - lambda .* (gradqx*u1_loc + gradqy*u2_loc);
            
            due1h_loc  = phiq*du1_loc;
            due2h_loc  = phiq*du2_loc;
            
            
            % Exact and approximated solutions at quadrature points
            ue1ex_loc = Data.ue_ex{1}(xq,yq)*Data.ue_t_ex{1}(time);
            ue2ex_loc = Data.ue_ex{2}(xq,yq)*Data.ue_t_ex{1}(time);
            due1ex_loc = Data.ue_ex{1}(xq,yq)*Data.due_t_ex{1}(time);
            due2ex_loc = Data.ue_ex{2}(xq,yq)*Data.due_t_ex{1}(time);
            
            % calcolare pressione esatta 
            
            
            % Fill the output Ue, dUe and Pe
            Ue(kinde + 1 : kinde + lqn,:) = [xq, yq, ue1h_loc, ue2h_loc, ue1ex_loc, ue2ex_loc];
            Pe(kinde + 1 : kinde + lqn,:) = [xq, yq, pe_loc];
            dUe(kinde + 1 : kinde + lqn,:) = [xq, yq, due1h_loc, due2h_loc, due1ex_loc, due2ex_loc];
            kinde = kinde + lqn ;
            
        end
    end
    
    
    
end

%% Output
Displacement.Ue  = Ue;
Displacement.UeTag = 'u_e';

Pressure.Ue = Pe;

Velocity.Ue  = dUe;
Velocity.UeTag = 'u_{e,t}';


