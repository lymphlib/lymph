%> @file  Evalshape2DTriangles.m
%> @author Mattia Corti, Paola F. Antonietti, Caterina Leimer Saglio
%> @date 19 October 2025 
%> @brief Construction of the the structures for the Jacobi basis functions for
%> triangular elements.
%> 
%> The function constructs the Jacobi basis functions \f$\varphi_j(x)\f$ and the 
%> corresponding gradients \f$\nabla\varphi_j(x)\f$, evaluated at the nodes of quadrature for
%> triangular elements.
%> 
%======================================================================
%> @section classEvalshape2DTriangles Class description
%======================================================================
%> @brief  Construction of the the structures for the Jacobi basis functions for
%> triangular elements.
%>
%> @param femregion     Structure containing all the information of the
%> about the finite element approximation. 
%> @param ElemIdx       Identifier of the mesh element.
%> @param Nodes         Coordinates of the nodes of quadrature.
%>
%> @retval dphiq        Basis functions evaluated at the quadrature nodes.
%> @retval Grad         Gradient of basis functions evaluated at the quadrature nodes.
%======================================================================
 
function [dphiq, Grad]= Evalshape2DTriangles(femregion, ElemIdx, Nodes)

    % Extraction of the degree of discretization
    N = femregion.degree(ElemIdx);

    loc_coord = femregion.coord(femregion.connectivity{ElemIdx},:);

    % Jacobian of elemental map
    BJ = [loc_coord(2,1)-loc_coord(1,1), loc_coord(3,1)-loc_coord(1,1); ...
          loc_coord(2,2)-loc_coord(1,2), loc_coord(3,2)-loc_coord(1,2)];

    % Translation vector
    Trans = [loc_coord(1,1); loc_coord(1,2)];

    % Inverse of the transformations
    BJ_inv    = [BJ(2,2) -BJ(1,2); -BJ(2,1) BJ(1,1)]/det(BJ);
    Trans_inv = [-BJ(2,2)*Trans(1)+BJ(1,2)*Trans(2);BJ(2,1)*Trans(1)-BJ(1,1)*Trans(2)]/det(BJ);

    % Computation of reference nodes
    NodesPhys = (BJ_inv * Nodes' + Trans_inv)';

    r = NodesPhys(:,1);
    s = NodesPhys(:,2);

    % Transfer of triangle coordinate from (r,s) to (a,b)
    a = zeros(size(s));
    a(s ~= 1) = (2*r(s ~= 1)./(1-s(s ~= 1))-1);
    a(s == 1) = - s(s == 1);
    b = 2*s-1;

    % Preallocation of the gradients
    Grad  = zeros(length(a),2,(N+1)*(N+2)/2);

    % First indexes vector
    idx1 = reshape(tril((1:N+1).*ones(N+1,N+1)),[1 (N+1)^2]);
    idx1(idx1==0) = [];
    idx1 = idx1 - 1;

    % Second indexes vector 
    idx2 = reshape(flip(triu((N+1:-1:1).*ones(N+1,N+1)),2)',[1 (N+1)^2]);
    idx2(idx2==0) = [];
    idx2 = idx2 - 1;

    if nargout == 1
            dphiq = FullSimplex2DP(a,b,idx1,idx2);
    else
            [dphiq, Grad(:,1,:), Grad(:,2,:)] = FullSimplex2DP(a,b,idx1,idx2);

            Grad = pagemtimes(Grad,BJ_inv);
    end
    
end
