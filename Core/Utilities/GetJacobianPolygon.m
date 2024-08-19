%> @file  GetJacobianPolygon.m
%> @author Mattia Corti
%> @date 11 November 2023 
%> @brief Jacobian map from reference to polygon vertices
%> 
%======================================================================
%> @section classDecr Class description
%======================================================================
%> @brief Jacobian map from reference to polygon vertices
%> 
%> The function constructs and evaluates the Jacobian of the
%> transformation in the polygon vertices.
%>
%>
%> @param LocCoord   Matrix containing the local coordinates listed as follows: 
%> \f[\begin{bmatrix}
%>        x_0 & y_0 \\
%>        x_1 & y_1 \\
%>        x_2 & y_2 \\
%>     \end{bmatrix} 
%>  \f]
%> @param v         Initial polygon vertices.
%>
%> @retval v        Polygon coordinates in reference configuration.
%> @retval Jdet     Determinant of the Jacobianof the transformation.
%> @retval BJ       Jacobian of the transformation.
%> @retval Trans    Translation vector of the transformation.


%======================================================================

function [v, Jdet, BJ, Trans] = GetJacobianPolygon(LocCoord, v)

        % Translation vector for the polygon vertices
        Trans = [LocCoord(1) + LocCoord(2);  LocCoord(3) + LocCoord(4)]/2;
        
        % Jacobian of the transformation
        BJ = 0.5 * [LocCoord(2) - LocCoord(1), 0; ...
                    0, LocCoord(4) - LocCoord(3)];
        
        % Determinant of the Jacobian of the transformation
        Jdet = det(BJ);

        % Polygon in the reference configuration
        v = BJ\(v - Trans);
end