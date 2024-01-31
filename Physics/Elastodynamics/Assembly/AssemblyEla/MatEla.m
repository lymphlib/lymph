%> @file  MatEla.m
%> @author Ilario Mazzieri
%> @date 16 April 2023
%> @brief Assembly of the matrices for the elastic problem \cite AM2017
%>
%==========================================================================
%> @section classMatEla Class description
%==========================================================================
%> @brief Assembly of the matrices for the elastic problem \cite AM2017
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Matrices.Ela = Elastic matrices (Mass, Stiffness,
%dG, projection matrix, etc)
%>
%==========================================================================

function [Matrices] = MatEla(Data, neighbor, femregion)
%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices
A = AllocateMatrixElaGlobal(femregion, neighbor);

%% Loop over the elements

% Visualization of computational progress
index_shift = 0;
id_shift = 0;

prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);

for ie = 1 : femregion.nel 

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
    [El] = AllocateMatrixElaLocalEl(femregion.nbases);
    
    % Check if the element is poroelastic
    if tag_ie == 'E'
        
        index_e = index - (femregion.ndof_p + femregion.ndof_a);

        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Evaluation of physical parameters
            par.rho_e = Data.rho_el{id_ie-id_shift}(xq,yq);
            par.vs    = Data.vs_el{id_ie-id_shift}(xq,yq);
            par.vp    = Data.vp_el{id_ie-id_shift}(xq,yq);
            par.zeta  = Data.zeta{id_ie-id_shift}(xq,yq);
            par.mu    = Data.mu_el{id_ie-id_shift}(xq,yq);
            par.lam   = Data.lam_el{id_ie-id_shift}(xq,yq);
            par.m     = 0;
            par.beta  = 0;
            
            % Construction of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
                    
            %% Local matrix assembling
            [El] = MakeLocMatrixElaElem(El, dx, phiq, gradqx, gradqy, par);
            
        end

        % Local matrices to global matrices
        [A] = LocalToGlobalMatrixElaElem(A, El, index_e);
        
        %% Boundary integrals and stabilization terms
        
        % Local matrix allocation for edges
        % [Edge] = AllocateMatrixElaLocalEdge(femregion.nbases, neighbor.nedges(ie));
        [Edge] = AllocateMatrixElaLocalEdge(femregion.nbases, neigh_ie_unq);

        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

        for iedg = 1 : neighbor.nedges(ie) % loop over faces
            
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
            
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            
            % Scaled weights
            ds = meshsize(iedg) * w_1D;
            
            % Extraction of normals to the face
            nx = normals(1,iedg);
            ny = normals(2,iedg);
           
            % Evaluation of physical parameters
            par.rho_e = Data.rho_el{id_ie-id_shift}(xq,yq);
            par.vs    = Data.vs_el{id_ie-id_shift}(xq,yq);
            par.vp    = Data.vp_el{id_ie-id_shift}(xq,yq);
            par.zeta  = Data.zeta{id_ie-id_shift}(xq,yq);
            par.mu    = Data.mu_el{id_ie-id_shift}(xq,yq);
            par.lam   = Data.lam_el{id_ie-id_shift}(xq,yq);
            par.m     = 0;
            par.beta  = 0;

            if strcmp(tag_neigh,'E') %elastic neighbor
                par.mu_n   = Data.mu_el{id_neigh-id_shift}(xq,yq);
                par.lam_n  = Data.lam_el{id_neigh-id_shift}(xq,yq);
                par.m_n    = 0;
                par.beta_n = 0;
            elseif strcmp(tag_neigh,'P') %poroelastic neighbor
                par.mu_n   = Data.mu{id_neigh}(xq,yq);
                par.lam_n  = Data.lam{id_neigh}(xq,yq);
                par.m_n    = Data.m{id_neigh}(xq,yq);
                par.beta_n = Data.beta{id_neigh}(xq,yq);
            elseif strcmp(tag_neigh,'NaN') %boundary edge
                par.mu_n   = Data.mu_el{id_ie-id_shift}(xq,yq);
                par.lam_n  = Data.lam_el{id_ie-id_shift}(xq,yq);
                par.m_n    = 0;
                par.beta_n = 0;
                %PhD thesis I.Mazzieri page 40
                par.c11 = (-par.mu ./ par.vs * ny * ny - (par.lam + 2*par.mu) ./ par.vp * nx * nx);
                par.c12 = ( par.mu ./ par.vs * nx * ny - (par.lam + 2*par.mu) ./ par.vp * nx * ny);
                par.c21 = ( par.mu ./ par.vs * nx * ny - (par.lam + 2*par.mu) ./ par.vp * nx * ny);
                par.c22 = (-par.mu ./ par.vs * nx * nx - (par.lam + 2*par.mu) ./ par.vp * ny * ny);               
                par.c3_11 =  (par.mu .* (2*par.vs-par.vp) ./ par.vs + (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp ) *nx*ny^2;
                par.c4_11 = -(par.mu .* (2*par.vs-par.vp) ./ par.vs + (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp ) *ny*nx^2;               
                par.c3_12 =  (par.mu .* (2*par.vs-par.vp) ./ par.vs * ny^3    - (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp *ny*nx^2 );
                par.c4_12 = -(par.mu .* (2*par.vs-par.vp) ./ par.vs * nx*ny^2 + (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp *nx^3 );                
                par.c3_21 = -(par.mu .* (2*par.vs-par.vp) ./ par.vs * ny*nx^2 + (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp *ny^3 );
                par.c4_21 =  (par.mu .* (2*par.vs-par.vp) ./ par.vs * nx^3    - (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp *nx*ny^2 );                
                par.c3_22 = -(par.mu .* (2*par.vs-par.vp) ./ par.vs + (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp ) *nx*ny^2;
                par.c4_22 =  (par.mu .* (2*par.vs-par.vp) ./ par.vs + (par.lam .* par.vs + 2*par.mu .* (par.vs-par.vp)) ./ par.vp ) *ny*nx^2;
            end

            par.lambda_ave = 2*par.lam .* par.lam_n ./ (par.lam + par.lam_n);
            par.mu_ave     = 2*par.mu .* par.mu_n ./ (par.mu + par.mu_n);
            par.harm_ave   = (par.lambda_ave + 2*par.mu_ave);
            par.m_ave      = 2*par.m .* par.m_n ./ (par.m + par.m_n);
            par.beta_ave   = 2*par.beta .* par.beta_n ./ (par.beta + par.beta_n);           
                        
            % Construction of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);            
            
            % Dirichlet boundary faces
            if neigh_ie(iedg) == -1
                
                [Edge] = MakeLocMatrixElaEdgeDirichlet(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, par, nx, ny, penalty_geom(iedg));
                
            % Absorbing boundary faces
            elseif neigh_ie(iedg) == -3 

                [Edge] = MakeLocMatrixElaEdgeAbsorbing(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, par);
                                             
            % Elastic neighbor
            elseif neigh_ie(iedg) >0 && femregion.tag(neigh_ie(iedg)) == 'E'
                
                % Element itself
                [Edge] = MakeLocMatrixElaEdgeEla(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, par, nx, ny, penalty_geom(iedg));
                
                % Construction of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);
                
                % Neighboring element
%                 [Edge] = MakeLocMatrixElaEdgeElaNeigh(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, phiedgeqneigh, par, nx, ny, penalty_geom(iedg), iedg);
                [Edge] = MakeLocMatrixElaEdgeElaNeigh(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, phiedgeqneigh, par, nx, ny, penalty_geom(iedg), idneigh);

              
            end
            
        % Local matrix to global matrix
        [A] = LocalToGlobalMatrixElaEdge(A, Edge, index_e);
            
        end

        % Local matrix to global matrix neighbor 
        % [A] = LocalToGlobalMatrixElaEdgeNeigh(A, Edge, index_e, neigh_ie, femregion.nbases, 0, femregion.nel_e);
        [A] = LocalToGlobalMatrixElaEdgeNeigh(A, Edge, index_e, neigh_ie_unq, femregion.nbases, 0, femregion.nel_e);

        
    end
    
end

%% Symmetric Interior Penalty Method

% Mass matrix elastic
M_P_rho   = [A.M1_P_rho, 0*A.M1_P_rho; 0*A.M1_P_rho, A.M1_P_rho];
% Projection matrix
MprjP   = [A.MPrjP_1, 0*A.MPrjP_1; 0*A.MPrjP_1, A.MPrjP_1];
% Damping matrix velocity
D = [A.D1, 0*A.D1; 0*A.D1, A.D1];
C = [A.C1, 0*A.C1; 0*A.C1, A.C1];

% Absorbing matrices
ABC_S = -[A.ABC_S1, A.ABC_S2; A.ABC_S3, A.ABC_S4];
ABC_R = -[A.ABC_R1, A.ABC_R2; A.ABC_R3, A.ABC_R4];

% Elastic dg matrix
V    = [A.V1,    A.V2;    A.V3,    A.V4];
IT_P = [A.IT1_P, A.IT2_P; A.IT3_P, A.IT4_P];
S_P  = [A.S1_P,  A.S2_P;  A.S3_P,  A.S4_P];


Matrices.Ela = struct( ...
    'A_E', V + S_P - IT_P - transpose(IT_P), ...
    'Dvel', D, ...
    'Ddis', C, ...
    'Svel', ABC_S, ...
    'Rdis', ABC_R, ...
    'MPrjP',MprjP,...
    'DGe', V + S_P,...
    'M_P_rho', M_P_rho);

fprintf('\n');

