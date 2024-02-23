%> @file  GetSolutionQuadPoints.m
%> @author Ilario Mazzieri
%> @date 16 April 2023
%> @brief Compute the solution at quadrature nodes.
%>
%==========================================================================
%> @section classGetSolutionQuadPoints Class description
%==========================================================================
%> @brief  Compute the solution at quadrature nodes.
%>
%> @param Data        Struct with problem's data
%> @param femregion   Struct for finite elements data
%> @param neighbor    Struct for neighboring elements data
%> @param Solutions   Struct with modal solutions
%>
%> @retval U   Solution array having the following structure
%>             U = [ x | y | uh(x,y) | uex(x,y)];
%>             (x,y) is the list of quadrature points
%>
%==========================================================================
function [U] = GetSolutionQuadPoints(Data,femregion,neighbor,Solutions)
%% Setup
% U = [ x | y | uh(x,y) | uex(x,y)];
U  = zeros(femregion.nqn.^2*sum(neighbor.nedges-2),4);
% U  = zeros(1,4);
i_shift = 0;

%% Quadrature values
[~, ~, ref_qNodes_2D, ~] = Quadrature(femregion.nqn);

%% Loop over the elements
% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie=1:femregion.nel
        
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
        % Local solution
        u_loc = Solutions.U(index);
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [~, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            lqn = length(xq);
         
            % Construction of the basis functions
            [phiq, ~, ~] = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Exact and approximated solutions at quadrature points
            uex_loc = Data.u_ex{1}(xq,yq);
            uh_loc  = phiq*u_loc;
            
            % Fill the output U
            U(i_shift + 1 : i_shift + lqn,:) = [xq, yq, uh_loc, uex_loc];            
            i_shift = i_shift + lqn;
            
        end
    end
    
    
    
end
