%> @file  ForcingHeat.m
%> @author Mattia Corti
%> @date 7 July 2023
%> @brief Assembly of the right-hand-side for the heat equation
%
%==========================================================================
%> @section classForcingHeat Class description
%==========================================================================
%> @brief            Assembly the RHS and boundary conditions for the heat equation
%>
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param t          Time associated to the forcing term
%>
%> @retval F         Right-hand-side vector
%>
%==========================================================================

function [F] = ForcingHeat(Data, neighbor, femregion, t)

    %% Quadrature values

    [ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

    %% Initialization of the forcing term

    F = zeros(femregion.ndof,1);


    %% Loop over the elements
    for ie = 1:femregion.nel

        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

        % Extraction of neighbor element and their edges
        neigh_ie      = neighbor.neigh{ie};
        neighedges_ie = neighbor.neighedges{ie};

        % Extraction of element geometrical information
        coords_ie          = femregion.coords_element{ie};
        [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);

        % Creation of the subtriangulation of the element
        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

        % Local forcing vector definition
        F_loc = zeros(femregion.nbases,1);

        if ~Data.homog_source_f

            for iTria = 1:size(Tria,1)

                 % Construction of Jacobian and quadrature nodes
                [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

                xq  = qNodes_2D(:,1);
                yq  = qNodes_2D(:,2);

                % Scaled weights
                dx = det(BJ) * w_2D;

                % Construction of the basis functions
                phiq = Evalshape2D(femregion, ie, qNodes_2D);

                % Evaluation of physical parameters
                mu   	= Data.mu{1}(xq,yq);
                f       = Data.source_f{1}(xq,yq,t);

                %% Vector assembling

                F_loc = F_loc + (dx.*phiq)'*f;

             end
        end

        if Data.TagApplyBCs == 1

            % Computation of all the penalty coefficients for the element ie
            [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

            % Loop over faces
            for iedg = 1 : neighbor.nedges(ie)

                % Extraction of the edge coordinates
                if iedg == neighbor.nedges(ie)
                    p1 = coords_ie(iedg,:);
                    p2 = coords_ie(1,:);
                else
                    p1 = coords_ie(iedg,:);
                    p2 = coords_ie(iedg+1,:);
                end

                % Construction of quadrature nodes on the face
                [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);

                xq = qNodes_1D(:,1);
                yq = qNodes_1D(:,2);

                % Scaled weights
                ds = meshsize(iedg) * w_1D;

                % Construction of the basis functions
                [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

                % Extraction of normals to the face
                nx = normals(1,iedg);
                ny = normals(2,iedg);

                % Dirichlet boundary faces
                if neigh_ie(iedg) == -1

                    % Physical parameters
                    mu     = Data.mu{1}(xq,yq);
                    gD     = Data.DirBC(xq,yq,t);

                    %% Vector assembling
                    F_loc = F_loc - (ds .* mu .* ( nx * gradedgeqx + ny * gradedgeqy))' * gD;
                    F_loc = F_loc + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * gD;

                end
            end
         end

        % Local vector to global vector
        F(index) = F_loc;
    end
end
