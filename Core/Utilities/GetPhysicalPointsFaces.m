
%> @file  GetPhysicalPointsFaces.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief Map from reference to physical points (face)
%>
%======================================================================
%> @section classGetPhysicalPointsFaces Class description
%======================================================================
%> @brief Map from reference to physical points (face)
%> 
%> The function for each face computes the Jacobian of the mapping 
%> from the reference points in \f$[-1,1]\f$ to the physical points.
%>
%>
%> @param LocCoord   Matrix containing the local coordinates listed as follows: 
%> \f[\begin{bmatrix}
%>        x_0 & y_0 \\
%>        x_1 & y_1 \\
%>     \end{bmatrix} 
%>  \f]
%> @param Nodes1D    Matrix containing the reference points for the evaluation.
%>
%> @retval PPhys1D   Evaluation of the physical points of the face.

%======================================================================
% for each face compute the mapping from the reference points to the
% physical points.

function [PPhys1D] = GetPhysicalPointsFaces(LocCoord, Nodes1D)

    % Jacobian of the transformation
    BJ_face = [LocCoord(2,1) - LocCoord(1,1), 0.5*LocCoord(2,1) - 0.5*LocCoord(1,1); ...
               LocCoord(2,2) - LocCoord(1,2), 0.5*LocCoord(2,2) - 0.5*LocCoord(1,2)];

    % Translation vector
    Trans = [LocCoord(1,1); LocCoord(1,2)];

    % Recover the dimensionality of the space (2)
    NodesFace = [Nodes1D zeros(length(Nodes1D), 1)]';

    % Physical points
    PPhys1D = (BJ_face * NodesFace + Trans)';

end


