%> @file  GetSolutionQuadPoints.m
%> @author Mattia Corti
%> @date 3 October 2023
%> @brief Evaluate the numerical and exact solution in the quadrature nodes
%>
%==========================================================================
%> @section classHeatGetSolutionQuadPoints Class description
%==========================================================================
%> @brief            Evaluate the exact solution constructing the vector
%> associated to the function in the modal basis
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param u_h        Numerical solution to evaluate in quadrature nodes
%> @param t          Time associated to the forcing term
%>
%> @retval Gc        Matrix containing in the first two columns the
%> quadrature nodes, in the third one the numerical solution and in the
%> fourth one the exact solution (if the last column exists):
%> Gc = [ x | y | u_h(x,y) | u_ex(x,y)]
%>
%==========================================================================

function [Gc] = GetSolutionQuadPoints(Data, femregion, neighbor, u_h, t)
    
    %% Quadrature values
    
    [~, ~, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

    %% Initialization of the forcing term
    Gc = zeros(sum(neighbor.nedges-2)*length(w_2D),4);
    
    %% Loop over the elements
    
    ContID = 0;

    for ie = 1:femregion.nel
        
        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        
        % Extraction of element geometrical information
        coords_ie          = femregion.coords_element{ie};

        % Creation of the subtriangulation of the element
        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
    
        for iTria = 1:size(Tria,1)
            
             % Construction of Jacobian and quadrature nodes
            [~, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
        
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Construction of the basis functions
            phiq = Evalshape2D(femregion, ie, qNodes_2D);
            
            %% Exact solution
                
            local_approx_u = phiq*u_h(index);
        
            if Data.PlotExact
                % Evaluation of exact solution
                local_exact_u = Data.u_ex(xq,yq,t);

                % Construction of the output structure
                Gc(ContID+1:ContID+length(xq),1:4) = [xq, yq, local_approx_u, local_exact_u];
            else

                % Construction of the output structure
                Gc(ContID+1:ContID+length(xq),1:3) = [xq, yq, local_approx_u];
            end

            ContID = ContID + length(xq);
        end

    end

end
