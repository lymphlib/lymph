%> @file  MatFluid.m
%> @author Ilario Mazzieri
%> @date 13 Febraury 2024
%> @brief Assembly of the matrices for the Stokes problem in pseudo-stress
%formulation
%>
%==========================================================================
%> @section classMatFluid Class description
%==========================================================================
%> @brief Assembly of the matrices for the Stokes problem in pseudo-stress
%formulation
%>
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%>
%> @retval Matrices  Matrices.Fluid = Fluid matrices (Mass, Stiffness,
%dG, projection matrix, etc)
%>
%==========================================================================

function [Matrices] = MatFluid(Data, neighbor, femregion)

%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices
A = AllocateMatrixFluidGlobal(femregion, neighbor);

%% Loop over the elements

% Visualization of computational progress
index_shift = 0;
id_shift = 0;

prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);

for ie = 1 : femregion.nel % loop over elements
    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);

    % Id and tag selection for elements
    id_ie  = femregion.id(ie);
    tag_ie = femregion.tag(ie);

    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    index_element = index_shift + (1:1:femregion.nedges(ie))';
    index_shift   = index_element(end);

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
    [El] = AllocateMatrixFluidLocalEl(femregion.nbases);

    % Check if the element is fluid
    if tag_ie == 'F'

        index_f = index;

        for iTria = 1:size(Tria,1)

            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Scaled weights
            dx = det(BJ) * w_2D;

            % Evaluation of physical parameters
            par.mu = Data.mu_f{id_ie-id_shift}(xq,yq);

            % Construction of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);


            %% Local matrix assembling
            [El] = MakeLocMatrixFluidElem(El, dx, phiq, gradqx, gradqy, par);

        end

        % Local matrices to global matrices
        [A] = LocalToGlobalMatrixFluidElem(A, El, index_f);

        %% Boundary integrals and stabilization terms

        % Local matrix allocation for edges
        % [Edge] = AllocateMatrixFluidLocalEdge(femregion.nbases, neighbor.nedges(ie));
        [Edge] = AllocateMatrixFluidLocalEdge(femregion.nbases, neigh_ie_unq);

        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);


        for iedg = 1 : neighbor.nedges(ie)

            % Extraction of tag and id of neighbor el
            id_el_neigh = neigh_ie(iedg);
            idneigh = (neigh_ie_unq == neighbor.neigh{ie}(iedg));

            if id_el_neigh > 0
                id_neigh  = femregion.id(id_el_neigh);
                tag_neigh = femregion.tag(id_el_neigh);
            else
                id_neigh  = 0;
                tag_neigh = 'NaN';
            end

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

            % Scaled weights
            ds = meshsize(iedg) * w_1D;

            % Extraction of normals to the face
            nx = normals(1,iedg);
            ny = normals(2,iedg);

            % Construction of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);


            %% Matrix assembling

            % Dirichlet boundary faces sigma*n = gn
            if neigh_ie(iedg) == -2

                [Edge] = MakeLocMatrixFluidEdgeDirichlet(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, nx, ny, penalty_geom(iedg));

                % Fluid neighbor
            elseif neigh_ie(iedg) >0 && femregion.tag(neigh_ie(iedg)) == 'F'

                % Element itself
                [Edge] = MakeLocMatrixFluidEdgeFluid(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, nx, ny, penalty_geom(iedg));

                % Construction of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);

                % Neighboring element
%                 [Edge] = MakeLocMatrixFluidEdgeFluidNeigh(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, phiedgeqneigh, nx, ny, penalty_geom(iedg), iedg);
                [Edge] = MakeLocMatrixFluidEdgeFluidNeigh(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, phiedgeqneigh, nx, ny, penalty_geom(iedg), idneigh);

            end

            % Local matrix to global matrix
            [A] = LocalToGlobalMatrixFluidEdge(A, Edge, index_f);

        end

        % Local matrix to global matrix neighbor
        % [A] = LocalToGlobalMatrixFluidEdgeNeigh(A, Edge, index_f, neigh_ie, femregion.nbases, femregion.nel_p, femregion.nel_p + femregion.nel_f);
        [A] = LocalToGlobalMatrixFluidEdgeNeigh(A, Edge, index_f, neigh_ie_unq, femregion.nbases, 0, femregion.nel_f);


    end

end

%% Symmetric Interior Penalty Method

% Mass matrices
MM_F   = [A.MF_2, A.Z,     A.Z,     A.Z;
    A.Z,     A.MF_2, A.Z,     A.Z;
    A.Z,     A.Z,     A.MF_2, A.Z;
    A.Z,     A.Z,     A.Z,     A.MF_2];

MM_DEV = [ A.Z,     A.Z,     A.Z,    -A.MF_2;
    A.Z,     A.MF_2,  A.Z,     A.Z;
    A.Z,     A.Z,     A.MF_2,  A.Z;
    -A.MF_2,  A.Z,     A.Z,     A.Z];

M_F = 0.5 * MM_F + 0.5 * MM_DEV;

% Projection matrix
MprjP   = [A.Mprj, A.Z,    A.Z,    A.Z;
    A.Z,    A.Mprj, A.Z,    A.Z;
    A.Z,    A.Z,    A.Mprj, A.Z;
    A.Z,    A.Z,    A.Z,    A.Mprj];

% Fluid dg matrix
BB   = [A.B1,   A.B2;   A.B3,   A.B4];
BBT  = [A.BT1,  A.BT2;  A.BT3,  A.BT4];
S_BB = [A.S1_B, A.S2_B; A.S3_B, A.S4_B];

ZZ = 0.*BB;

B   = [BB   ZZ; ZZ BB  ];
BT  = [BBT  ZZ; ZZ BBT ];
S_B = [S_BB ZZ; ZZ S_BB];


Matrices.Fluid = struct(...
    'MPrj', A.Mprj,...
    'A_F', B - BT - transpose(BT) + S_B,...
    'MPrjP',MprjP,...
    'DGf',B + S_B,...
    'M_F', M_F);

fprintf('\n');


