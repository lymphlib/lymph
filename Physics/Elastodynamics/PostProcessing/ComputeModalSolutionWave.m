%> @file  ComputeModalSolutionWave.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Compute modal expansions for elastic displacement
%>
%==========================================================================
%> @section classComputeModalSolutionWave Class description
%> @brief  Compute Modal Solution for the Wave problem
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%
%> @retval ue       Vector with modal expansions
%>
%==========================================================================

function [ue] = ComputeModalSolutionWave(Data,femregion)
%% Quadrature values
[~, ~, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization
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
    if tag_ie == 'E'
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        % Inizialize local silution vector
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

ue = [ue1; ue2];

fprintf('\n');

