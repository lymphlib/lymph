%> @file  Evalshape2D.m
%> @author Mattia Corti, Paola F. Antonietti, Caterina Leimer Saglio
%> @date 28 October 2025
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

function [phiq, gradqx, gradqy, lapqxx, lapqxy, lapqyx, lapqyy] = Evalshape2D(femregion, ElemIdx, Nodes)
     
    % Shape of mesh handling
    if all(femregion.nedges == 3)

        % Triangular mesh
        Evalshape2DElem = @Evalshape2DTriangles;
    else

        % Polygonal mesh
        Evalshape2DElem = @Evalshape2DPolygons;
    end

    switch nargout 
        
        % Computation of bases
        case 1
            phiq = Evalshape2DElem(femregion, ElemIdx, Nodes);
     
        % Computation of bases and gradients 
        case 3

            [phiq, Grad] = Evalshape2DElem(femregion, ElemIdx, Nodes);
       
            % Gradient reshape
            [gradqx, gradqy] = extractGrad(Grad, Nodes);
        
        % Computation of bases, gradients, and laplacians
        otherwise

            if all(femregion.nedges == 3)
                error("Laplacian computation not implemented for simplexes")
            end
            
            [phiq, Grad, Lap] = Evalshape2DElem(femregion, ElemIdx, Nodes);

            % Gradient reshape
            [gradqx, gradqy] = extractGrad(Grad, Nodes);

            % Laplacian reshape
            [lapqxx, lapqxy, lapqyx, lapqyy] = extractLap(Lap, Nodes);

    end

end



%% Internal function for gradient extraction
function [gradqx, gradqy] = extractGrad(Grad, Nodes)

    gradqx = squeeze(Grad(:,1,:));
    gradqy = squeeze(Grad(:,2,:));

    if size(Nodes,1) == 1
        gradqx = gradqx';
        gradqy = gradqy';
    end
end



%% Internal function for laplacian extraction
function [lapqxx, lapqxy, lapqyx, lapqyy] = extractLap(Lap, Nodes)

    lapqxx = squeeze(Lap(:,1,:));
    lapqxy = squeeze(Lap(:,2,:));
    lapqyx = squeeze(Lap(:,3,:));
    lapqyy = squeeze(Lap(:,4,:));

    if size(Nodes,1) == 1
        lapqxx = lapqxx';
        lapqxy = lapqxy';
        lapqyx = lapqyx';
        lapqyy = lapqyy';
    end

end
