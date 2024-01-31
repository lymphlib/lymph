%> @file  EvaluateSolution.m
%> @author Mattia Corti
%> @date 3 October 2023
%> @brief Evaluate the exact solution constructing the vector associated
%>
%==========================================================================
%> @section classEvaluateSolution Class description
%==========================================================================
%> @brief            Evaluate the exact solution constructing the vector
%> associated to the function in the modal basis
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param t          Time associated to the forcing term
%>
%> @retval u         Exact-solution associated vector
%>
%==========================================================================

function [u] = EvaluateSolution(Data, femregion, t)

    %% Quadrature values
    
    [~, ~, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);
    

    %% Initialization of the forcing term

    u = zeros(femregion.ndof,1);

    
    %% Loop over the elements
       
    for ie = 1:femregion.nel
        
        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        
        % Extraction of element geometrical information
        coords_ie          = femregion.coords_element{ie};
        
        % Creation of the subtriangulation of the element
        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
    
        % Local forcing vector definition
        u_loc = zeros(femregion.nbases,1);

        for iTria = 1:size(Tria,1)
                 
             % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
        
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Scaled weights 
            dx = det(BJ) * w_2D;

            % Construction of the basis functions
            phiq = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Evaluation of physical parameters 
            u_ex   = Data.u_ex(xq,yq,t);
            
            %% Vector assembling

            u_loc = u_loc + (dx .* phiq)'*u_ex;

         end

         % Local vector to global vector
         u(index) = u_loc;
    end
    
end
