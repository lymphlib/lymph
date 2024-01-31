%> @file  GetJacobianPhysicalPoints.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief Jacobian map from reference to physical points 
%> 
%======================================================================
%> @section classGetJacobianPhysicalPoints Class description
%======================================================================
%> @brief Jacobian map from reference to physical points 
%> 
%> The function constructs and evaluates the Jacobian of the
%> transformation in the physical points.
%>
%>
%> @param LocCoord   Matrix containing the local coordinates listed as follows: 
%> \f[\begin{bmatrix}
%>        x_0 & y_0 \\
%>        x_1 & y_1 \\
%>        x_2 & y_2 \\
%>     \end{bmatrix} 
%>  \f]
%> @param Nodes2D    Matrix containing the reference points for the evaluation.
%>
%> @retval BJ         Jacobian matrix:
%> \f[\begin{bmatrix}
%>        x_1 - x_0 & x_2 - x_0 \\
%>        y_1 - y_0 & y_2 - y_0 \\
%>     \end{bmatrix} 
%>  \f]
%> @retval PPhys2D   Evaluation of the jacobian at the physical points.

%======================================================================

function [BJ, PPhys2D] = GetJacobianPhysicalPoints(LocCoord, Nodes2D)

    % Jacobian of the Transformation
    BJ = [LocCoord(2,1) - LocCoord(1,1), LocCoord(3,1) - LocCoord(1,1); ...
          LocCoord(2,2) - LocCoord(1,2), LocCoord(3,2) - LocCoord(1,2)];
    
    % Translation vector
    Trans = [LocCoord(1,1); LocCoord(1,2)];
    
    % Physical points
    PPhys2D = (BJ * Nodes2D' + Trans)';
end