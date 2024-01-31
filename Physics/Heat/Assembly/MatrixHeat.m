%> @file  MatrixHeat.m
%> @author Mattia Corti
%> @date 3 October 2023
%> @brief Assembly of the matrices for the heat equation
%
%==========================================================================
%> @section classMatrixHeat Class description
%==========================================================================
%> @brief            Assembly of the mass and stiffness matrices for the heat equation.
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%>
%> @retval Matrices  Matrices.M       = Mass Matrix
%>                   Matrices.M_prj   = Mass Matrix (for projection)
%>                   Matrices.A       = Stiffness matrix
%>                   Matrices.dGA     = dG matrix for error analysis
%>
%==========================================================================

function [Matrices] = MatrixHeat(Data, neighbor, femregion)

    %% Quadrature values

    [ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

    %% Initialization of the matrices

    M     = spalloc(femregion.ndof,femregion.ndof,femregion.nel*femregion.nbases^2);      % Mass Matrix

    M_prj = spalloc(femregion.ndof,femregion.ndof,femregion.nel*femregion.nbases^2);      % Mass projection Matrix

    A     = spalloc(femregion.ndof,femregion.ndof,femregion.nel*femregion.nbases^2);      % Stiffness Matrix

    IA    = spalloc(femregion.ndof,femregion.ndof,femregion.nbases^2*sum(neighbor.nedges+1));     % DG-Stiffness Matrix

    SA    = spalloc(femregion.ndof,femregion.ndof,femregion.nbases^2*sum(neighbor.nedges+1));     % DG-Stiffness Matrix


    %% Loop over the elements

    % Visualization of computational progress
    prog = 0;
    fprintf(1,'Computation Progress: %3d%%\n',prog);

    for ie = 1:femregion.nel

        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);

        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

        % Extraction of neighbor element and their edges
        neigh_ie      = neighbor.neigh{ie};
        neigh_ie_unq  = unique(neighbor.neigh{ie});
        neighedges_ie = neighbor.neighedges{ie};

        % Extraction of element geometrical information
        coords_ie          = femregion.coords_element{ie};
        [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);

        % Creation of the subtriangulation of the element
        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

        % Local matrices definitions
        A_loc     = zeros(femregion.nbases,femregion.nbases);
        M_loc     = zeros(femregion.nbases,femregion.nbases);
        M_prj_loc = zeros(femregion.nbases,femregion.nbases);

        % Loop over the subtriangulation
        for iTria = 1:size(Tria,1)

            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Scaled weights
            dx = det(BJ) * w_2D;

            % Evaluation of physical parameters
            mu    = Data.mu{1}(xq,yq);
	        sigma = Data.sigma{1}(xq,yq);

            % Construction of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

            %% Matrix assembling
            M_prj_loc = M_prj_loc + (dx .* phiq)' * phiq;      % Mass Projection Matrix

            M_loc = M_loc + (dx .* sigma .* phiq)' * phiq;      % Mass Matrix

            A_loc = A_loc + (dx .* (mu .* gradqx))' * gradqx + (dx .* (mu .* gradqy))' * gradqy;    % Stiffness Matrix

        end

        A(index,index)        = A_loc;
        M(index,index)        = M_loc;
        M_prj(index,index)    = M_prj_loc;

        %% Boundary integrals and stabilization terms

        % Local matrix allocation
        SA_loc = zeros(femregion.nbases,femregion.nbases);
        IA_loc = zeros(femregion.nbases,femregion.nbases);

        % Local matrix allocation for neighbors
        IAN_loc = zeros(femregion.nbases, femregion.nbases, length(neigh_ie_unq));
        SAN_loc = zeros(femregion.nbases, femregion.nbases, length(neigh_ie_unq));

         % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

        % Loop over faces
        for iedg = 1 : neighbor.nedges(ie)

            % Extraction of the id of the neighboring element in the matrices IAN and SAN
            idneigh = (neigh_ie_unq == neighbor.neigh{ie}(iedg));

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

            % Evaluation of physical parameters
            mu    = Data.mu{1}(xq,yq);
            sigma = Data.sigma{1}(xq,yq);

            % Construction of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

            % Extraction of normals to the face
            nx = normals(1,iedg);
            ny = normals(2,iedg);

            %% Matrix assembling

            % Dirichlet boundary faces
            if neigh_ie(iedg) == -1

                IA_loc = IA_loc + (ds .* mu .* (nx * gradedgeqx + ny * gradedgeqy))' * phiedgeq;
                SA_loc = SA_loc + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeq;

            % Internal faces
            elseif neigh_ie(iedg) > 0

                % Element itself
                IA_loc = IA_loc + 0.5 * (ds .* mu .* (nx * gradedgeqx + ny * gradedgeqy))' * phiedgeq;
                SA_loc = SA_loc + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeq;

                % Construction of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);

                % Neighboring element

                % Extracellular component
                IAN_loc(:,:,idneigh) = IAN_loc(:,:,idneigh) - 0.5 * (ds .* mu .* ( nx * gradedgeqx +  ny * gradedgeqy))' * phiedgeqneigh;
                SAN_loc(:,:,idneigh) = SAN_loc(:,:,idneigh) - penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeqneigh;

            end

            % Local matrix to global matrix
            IA(index,index) = IA_loc;
            SA(index,index) = SA_loc;
        end

        % Assembling I-matrices
        [IA] = AssembleNeighEl(IA, index, neigh_ie_unq, IAN_loc, femregion.nbases, 1, femregion.nel);
        [SA] = AssembleNeighEl(SA, index, neigh_ie_unq, SAN_loc, femregion.nbases, 1, femregion.nel);

    end

    %% Symmetric Interior Penalty Method

    dGA = A + SA;                   % For Error Computation in DG Norm

    A = dGA - IA - transpose(IA);   % DG Stiffness Matrix

    % Build the Matrices structure
    Matrices = struct('M',          M, ...
                      'M_prj',      M_prj, ...
                      'A',          A, ...
                      'dGA',        dGA);

end
