%> @file  MatrixLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Assembly of the matrices for the Poisson's problem
%>
%==========================================================================
%> @section classMatrixLaplacian Class description
%==========================================================================
%> @brief            Assembly of the mass and stifness matrices for the
%Poisson's problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Matrices.Mprj = Mass Matrix;
%>                   Matrices.A    = Stiffness matrix;
%>                   Matrices.dGA  = dG matrix for error analysis;
%>
%==========================================================================

function [Matrices] = MatrixLaplacian(Data, neighbor, femregion)

%% Quadrature values

[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices

% \int_{\Omega}  ( u.v ) dx
Mprj = spalloc(femregion.ndof, femregion.ndof,femregion.nel*femregion.nbases^2);

% \int_{\Omega} ( D grad(u).grad(v) ) dx
A = spalloc(femregion.ndof, femregion.ndof,femregion.nel*femregion.nbases^2);

% \int_{E_h} ( {D grad(u)}.[v] ) ds
IA = spalloc(femregion.ndof,femregion.ndof,femregion.nbases^2*sum(neighbor.nedges+1));

% \int_{E_h} penalty_E ( h_s^(-1) [u].[v] ) ds
SA = spalloc(femregion.ndof,femregion.ndof,femregion.nbases^2*sum(neighbor.nedges+1));


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
    neighedges_ie = neighbor.neighedges{ie};
    
    % Extraction of element geometrical information
    coords_ie          = femregion.coords_element{ie};
    [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);
    
    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
    
    % Local matrices definitions
    Mprj_loc = zeros(femregion.nbases, femregion.nbases);
    A_loc    = zeros(femregion.nbases, femregion.nbases);
    
    
    % Loop over the subtriangulation
    for iTria = 1:size(Tria,1)
        
        % Construction of Jacobian and quadrature nodes
        [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
        
        xq  = qNodes_2D(:,1);
        yq  = qNodes_2D(:,2);
        
        % Scaled weights
        dx = det(BJ) * w_2D;
        
        % Evaluation of physical parameters
        mu = Data.mu{1}(xq,yq);
        
        % Construction of the basis functions
        [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
        
        %% Matrix assembling
        
        Mprj_loc = Mprj_loc + (dx.*phiq)'*phiq;
        A_loc = A_loc + (dx .* (mu .* gradqx))' * gradqx + (dx .* (mu .* gradqy))' * gradqy;
        
    end
    
    % Local matrices to global matrices
    Mprj(index,index) = Mprj_loc;
    A(index,index)    = A_loc;
    
    %% Boundary integrals and stabilization terms
    
    % Local matrix allocation
    IA_loc = zeros(femregion.nbases,femregion.nbases);
    SA_loc = zeros(femregion.nbases,femregion.nbases);
    
    % Local matrix allocation for neighbors
    IAN_loc = zeros(femregion.nbases, femregion.nbases, neighbor.nedges(ie));
    SAN_loc = zeros(femregion.nbases, femregion.nbases, neighbor.nedges(ie));
    
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
        
        % Evaluation of physical parameters
        mu = Data.mu{1}(xq,yq);
        
        % Construction of the basis functions
        [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);
        
        % Extraction of normals to the face
        nx = normals(1,iedg);
        ny = normals(2,iedg);
        
        %% Matrix assembling
        
        % Dirichlet boundary faces
        if neigh_ie(iedg) == -1
            
            % IA_ij = SUM_{q} ds(x_q) * mu(x_q) * grad(phi_i(x_q)).n(x_q) * phi_j(x_q)
            IA_loc = IA_loc + (ds .* mu .* ( nx * gradedgeqx + ny * gradedgeqy))' * phiedgeq;
            SA_loc = SA_loc + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeq;
            
            % Internal faces
        elseif neigh_ie(iedg) > 0
            
            % Element itself
            IA_loc = IA_loc + 0.5 * (ds .* mu .* ( nx * gradedgeqx + ny * gradedgeqy))' * phiedgeq;
            SA_loc = SA_loc + (ds .* (mu * penalty_geom(iedg)) .* phiedgeq)' * phiedgeq;
            
            % Construction of the basis functions for the neighbor
            phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);
            
            % Neighboring element
            IAN_loc(:,:,iedg) = IAN_loc(:,:,iedg) - 0.5 * (ds .* (mu .* ( nx * gradedgeqx +  ny * gradedgeqy)))' * phiedgeqneigh;
            SAN_loc(:,:,iedg) = SAN_loc(:,:,iedg) - (ds .* (mu * penalty_geom(iedg)) .* phiedgeq)' * phiedgeqneigh;
            
        end
        
        % Local matrix to global matrix
        IA(index,index) = IA_loc;
        SA(index,index) = SA_loc;
        
    end
    
    % Assembling boundary DG matrices (Interior penalty and stabilization)
    [IA] = AssembleNeighEl(IA, index, neigh_ie, IAN_loc, femregion.nbases, 1, femregion.nel);
    [SA] = AssembleNeighEl(SA, index, neigh_ie, SAN_loc, femregion.nbases, 1, femregion.nel);
    
    
end

%% Symmetric Interior Penalty Method

dGA = A + SA; % For error computation in dG-norm

A = dGA - IA - transpose(IA); % dG stiffness matrix

% Build the Matrices structure
Matrices = struct('Mprj', Mprj, ...
    'A',    A,    ...
    'dGA',  dGA);

end
