%> @file  Evalshape2DTriangles.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 26 March 2023 
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
    N = femregion.degree;

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
    %a = (s ~= 1).*(2*r./(1-s)-1) - (s == 1);
    a = zeros(size(s));
    a(s ~= 1) = (2*r(s ~= 1)./(1-s(s ~= 1))-1);
    a(s == 1) = - s(s == 1);
    b = 2*s-1;

    % Basis functions reconstruction
    sk = 1;

    dphiq = zeros(length(a),(N+1)*(N+2)/2);
    dpsi  = zeros(length(a),(N+1)*(N+2)/2,2);
    
    for i=0:N

        for j=0:N - i

            dphiq(:,sk) = Simplex2DP(a,b,i,j);
            [dpsi(:,sk,1), dpsi(:,sk,2)] = GradSimplex2DP(a,b,i,j);
            sk = sk + 1;

        end
        
    end

    % Gradient basis functions
    Grad(:,1,:) = dpsi(:,:,1);
    Grad(:,2,:) = dpsi(:,:,2);

    for j = 1:size(Grad,3)
        Grad(:,:,j) = Grad(:,:,j)*BJ_inv;
    end



