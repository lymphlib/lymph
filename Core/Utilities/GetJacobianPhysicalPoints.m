%> @file  GetJacobianPhysicalPoints.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 19 October 2025 
%> @brief Jacobian maps from reference to physical points for all the
%> subtriangles of the element ie
%> 
%======================================================================
%> @section classGetJacobianPhysicalPoints Class description
%======================================================================
%> @brief Jacobian map from reference to physical points for all the
%> subtriangles of the element ie 
%> 
%> The function constructs and evaluates the Jacobian of the
%> transformation in the physical points.
%>
%>
%> @param coords_ie  Matrix containing the coordinates of the considered
%> element.
%> @param Tria       Constructed Delaunay subtriangulation of the element.
%> @param Nodes2D    Matrix containing the reference points for the evaluation.
%>
%> @retval detBJ     Determinant of the Jacobian matrix:
%> \f[\begin{bmatrix}
%>        x_1 - x_0 & x_2 - x_0 \\
%>        y_1 - y_0 & y_2 - y_0 \\
%>     \end{bmatrix} 
%>  \f]
%> @retval PPhys2D   Evaluation of the jacobian at the physical points.

%======================================================================

function [detBJ, PPhys2D] = GetJacobianPhysicalPoints(coords_ie, Tria, Nodes2D)

    Loc_coords_x = reshape(coords_ie(reshape(Tria',1,size(Tria,1)*size(Tria,2)),1)',3,1,[]);
    Loc_coords_y = reshape(coords_ie(reshape(Tria',1,size(Tria,1)*size(Tria,2)),2)',3,1,[]);

    % Jacobian of the transformation
    BJ = [Loc_coords_x(2,:,:)-Loc_coords_x(1,:,:) Loc_coords_x(3,:,:)-Loc_coords_x(1,:,:); 
          Loc_coords_y(2,:,:)-Loc_coords_y(1,:,:) Loc_coords_y(3,:,:)-Loc_coords_y(1,:,:)];

    % Determinant of the transformation
    detBJ = reshape(BJ(1,1,:).*BJ(2,2,:)-BJ(2,1,:).*BJ(1,2,:),1,size(BJ,3));

    % Translation vector
    Trans = [Loc_coords_x(1,:,:); Loc_coords_y(1,:,:)];

    % Physical points
    PPhys2D = reshape((pagemtimes(BJ,Nodes2D') + Trans),2,[])';
    
end