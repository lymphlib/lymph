%> @file  Evalshape2D.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 7 October 2024
%> @brief Construction of the structures for the basis functions.
%> 
%> The function constructs the basis functions \f$\varphi_j(x)\f$ and the 
%> corresponding gradients \f$\nabla\varphi_j(x)\f$, evaluated at the nodes of quadrature,
%> by recalling the correct construction algorithm for the mesh elements.
%> 
%======================================================================
%> @section classEvalshape2D Class description
%======================================================================
%> @brief  Construction of the the structures for the basis functions.
%>
%> @param femregion     Structure containing all the information of the
%> about the finite element approximation. 
%> @param ElemIdx       Identifier of the mesh element.
%> @param Nodes         Coordinates of the nodes of quadrature.
%>
%> @retval phiq         Basis functions evaluated at the quadrature nodes.
%> @retval gradqx       Gradient (x-component) of basis functions evaluated at the quadrature nodes.
%> @retval gradqy       Gradient (y-component) of basis functions evaluated at the quadrature nodes.
%======================================================================

function [phiq, gradqx, gradqy] = Evalshape2D(femregion, ElemIdx, Nodes)
    
    % Control if gradients are needed or not

    if nargout == 1
        % Shape of mesh handling
        if all(femregion.nedges == 3)
            % Triangular mesh
            phiq = Evalshape2DTriangles(femregion, ElemIdx, Nodes);
        else
            % Polygonal mesh
            phiq = Evalshape2DPolygons(femregion, ElemIdx, Nodes);
        end
     
    else

        % Shape of mesh handling
        if all(femregion.nedges == 3)
            % Triangular mesh
            [phiq, Grad] = Evalshape2DTriangles(femregion, ElemIdx, Nodes);
        else
            % Polygonal mesh
            [phiq, Grad] = Evalshape2DPolygons(femregion, ElemIdx, Nodes);
        end
       
        gradqx = squeeze(Grad(:,1,:));
        gradqy = squeeze(Grad(:,2,:));

        if size(Nodes,1) == 1
            gradqx = gradqx';
            gradqy = gradqy';
        end
        
    end

end
