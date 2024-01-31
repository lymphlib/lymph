%> @file  Evalshape2DPolygons.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 25 February 2023 
%> @brief Construction of the the structures for the Legendre basis functions for
%> polygonal elements.
%> 
%> The function constructs the Legendre basis functions \f$\varphi_j(x)\f$ and the 
%> corresponding gradients \f$\nabla\varphi_j(x)\f$, evaluated at the nodes of quadrature for
%> polygonal elements.
%> 
%======================================================================
%> @section classEvalshape2DPolygons Class description
%======================================================================
%> @brief  Construction of the the structures for the Legendre basis functions for
%> polygonal elements.
%>
%> @param femregion     Structure containing all the information of the
%> about the finite element approximation. 
%> @param ElemIdx       Identifier of the mesh element.
%> @param Nodes         Coordinates of the nodes of quadrature.
%>
%> @retval dphiq        Basis functions evaluated at the quadrature nodes.
%> @retval Grad         Gradient of basis functions evaluated at the quadrature nodes.
%======================================================================

function [dphiq, Grad] = Evalshape2DPolygons(femregion, ElemIdx, Nodes)

    % Extraction of the degree of discretization
    N = femregion.degree;

    BBox = femregion.bbox(ElemIdx,:);

    % Jacobian of elemental map
    BJ = 0.5 * [BBox(2) - BBox(1), 0; 0, BBox(4) - BBox(3)];       

    % Translation vector
    Trans = 0.5 * [BBox(1) + BBox(2); BBox(3) + BBox(4)];   
    
    % Inverse of the transformations
    BJ_inv    = [BJ(2,2) -BJ(1,2); -BJ(2,1) BJ(1,1)]/det(BJ);
    Trans_inv = [-BJ(2,2)*Trans(1)+BJ(1,2)*Trans(2);BJ(2,1)*Trans(1)-BJ(1,1)*Trans(2)]/det(BJ);

    % Computation of reference nodes
    NodesPhys = (BJ_inv * Nodes' + Trans_inv)';

    a = NodesPhys(:,1);
    b = NodesPhys(:,2);
    
    % First indexes vector
    idx1 = reshape(tril((1:N+1).*ones(N+1,N+1)),[1 (N+1)^2]);
    idx1(idx1==0) = [];
    idx1 = idx1 - 1;

    % Second indexes vector 
    idx2 = reshape(flip(triu((N+1:-1:1).*ones(N+1,N+1)),2)',[1 (N+1)^2]);
    idx2(idx2==0) = [];
    idx2 = idx2 - 1;

    if nargout == 1

        % Basis functions reconstruction
        dphiq = FullHypercube2DP(a,b,idx1,idx2);

    else
       
        % Basis functions reconstruction
        [dphiq, Grad(:,1,:), Grad(:,2,:)] = FullHypercube2DP(a,b,idx1,idx2);
        
        for j = 1:size(Grad,3)
            Grad(:,:,j)= Grad(:,:,j)*BJ_inv;
        end

    end

end

