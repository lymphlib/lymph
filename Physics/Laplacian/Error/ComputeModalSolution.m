%> @file  ComputeModalSolution.m
%> @author Ilario Mazzieri
%> @date 16 April 2023
%> @brief Compute modal coefficient of the exact solution
%>
%==========================================================================
%> @section classComputeModalSolution Class description
%==========================================================================
%> @brief Compute modal coefficient of the exact solution
%
%> @param Data        Struct with problem's data
%> @param femregion   Structure containing all the information 
%> about the finite element approximation
%
%> @retval u_mod      Modal coefficients of the exact solution 
%>
%==========================================================================



function [u_mod] = ComputeModalSolution(Data,femregion)
%% Quadrature values
[~, ~, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization
u_mod = zeros(femregion.ndof,1);

%% Loop over the elements
% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie = 1:femregion.nel
    
    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    
    % Extraction of element geometrical information
    coords_ie = femregion.coords_element{ie};
    
    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
    
    
    if femregion.tag(ie) == 'L'
        % Inizialize local silution vector
        u_loc = 0;
        
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
            uex_loc = Data.u_ex{1}(xq,yq);
            u_loc = u_loc + (dx.*phiq)'*uex_loc;
            
        end
        u_mod(index) = u_loc;
        
    end
end

fprintf('\n');

