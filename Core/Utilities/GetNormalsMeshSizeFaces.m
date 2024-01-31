%> @file  GetNormalsMeshSizeFaces.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 26 March 2023 
%> @brief Computation of face diameter and outward normal
%> 
%======================================================================
%> @section classGetNormalsMeshSizeFaces Class description
%======================================================================
%> @brief Computation of face diameter and outward normal
%> 
%> This function computes, for each face of a given element, the diameter 
%> of the face and the corresponding unit outward normal
%>
%>
%> @param LocCoord   Matrix containing the local coordinates of the element vertices listed as follows: 
%> \f[\begin{bmatrix}
%>        x_0 & y_0 \\
%>        x_1 & y_1 \\
%>        \vdots & \vdots \\
%>        x_N & y_N \\
%>     \end{bmatrix} 
%>  \f]
%>
%> @retval Normal    Matrix containing the outward normals of each face of
%> the element.
%> 
%> @retval MeshSize  Vector containing the diameter of each element.

%======================================================================

function [Normal, MeshSize] = GetNormalsMeshSizeFaces(LocCoord)

    % Repeat the first element
    LocCoord = [LocCoord; LocCoord(1,:)];

    % Computing the unit normal for each face
    nn = [(LocCoord(2:end,2)-LocCoord(1:end-1,2))';(LocCoord(1:end-1,1)-LocCoord(2:end,1))'];
    Normal = nn./vecnorm(nn);
    
    % Computing the size of each face
    MeshSize = (vecnorm(nn))';

end
