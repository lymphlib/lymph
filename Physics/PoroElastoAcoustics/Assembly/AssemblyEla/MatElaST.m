%> @file  MatElaST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 6 August 2024
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
%> @retval Matrices  Elastic matrices (Mass, Stiffness, dG, projection matrix, etc)
%>
%==========================================================================

function [Matrices] = MatElaST(Data, neighbor, femregion)
%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices
max_nedges = max(neighbor.nedges+1);

ii_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
jj_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
ii_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));
jj_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));

%% Initialization of the matrices
El.V1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
El.V2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
El.V3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
El.V4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));

El.M1_P_rho_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
El.MPrjP_1_loc   = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
                  
El.D1_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
El.C1_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));

Edge.ABC_R1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
Edge.ABC_R2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
Edge.ABC_R3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
Edge.ABC_R4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));

Edge.ABC_S1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
Edge.ABC_S2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
Edge.ABC_S3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));
Edge.ABC_S4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_e),femregion.nbases,femregion.nbases,ones(1,femregion.nel_e));

Edge.S1_P_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));
Edge.S4_P_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));

Edge.IT1_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));
Edge.IT2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));
Edge.IT3_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));
Edge.IT4_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_e),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_e));
 
%% Loop over the elements

% Visualization of computational progress
index_shift = 0;
id_shift = max([Data.TagElAcu,Data.TagElPoro]);
if(isempty(id_shift)); id_shift = 0; end
nel_sh = (femregion.nel_a + femregion.nel_p);

prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);

for ie = 1 : femregion.nel

    ie_sh = ie - nel_sh;

    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    % Id and tag selection for elements
    id_ie  = femregion.id(ie);
    tag_ie = femregion.tag(ie);

    % Check if the element is poroelastic
    if tag_ie == 'E'

        % Selection of the matrix positions associated to element ie
        index = (ie_sh-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        index_element = index_shift + (1:1:femregion.nedges(ie))';
        index_shift   = index_element(end);
        ii_index{ie_sh} = repmat(index, 1, femregion.nbases);
        jj_index{ie_sh} = repmat(index', femregion.nbases, 1);

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
            El.V1_loc{ie_sh} = El.V1_loc{ie_sh} + (dx .* ((par.lam+2*par.mu) .* gradqx))' * gradqx + (dx .* (par.mu .* gradqy))' * gradqy;
            El.V2_loc{ie_sh} = El.V2_loc{ie_sh} + (dx .* (par.lam        .* gradqx))' * gradqy + (dx .* (par.mu .* gradqy))' * gradqx;
            El.V3_loc{ie_sh} = El.V3_loc{ie_sh} + (dx .* (par.lam        .* gradqy))' * gradqx + (dx .* (par.mu .* gradqx))' * gradqy;
            El.V4_loc{ie_sh} = El.V4_loc{ie_sh} + (dx .* ((par.lam+2*par.mu) .* gradqy))' * gradqy + (dx .* (par.mu .* gradqx))' * gradqx;


            El.M1_P_rho_loc{ie_sh} = El.M1_P_rho_loc{ie_sh}  + (dx .* (par.rho_e .* phiq))' * phiq;
            El.MPrjP_1_loc{ie_sh}  = El.MPrjP_1_loc{ie_sh}   + (dx .* phiq)' * phiq;

            El.D1_loc{ie_sh} = El.D1_loc{ie_sh}  + (dx .* (2 * par.rho_e .* par.zeta .* phiq))' * phiq;
            El.C1_loc{ie_sh} = El.C1_loc{ie_sh} + (dx .* (par.rho_e .* par.zeta.^2  .* phiq))' * phiq;
            
        end
 
        %% Boundary integrals and stabilization terms
        
        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

        for iedg = 1 : neighbor.nedges(ie) % loop over faces
            
            % Extraction of tag and id of neighbor el
            id_el_neigh = neigh_ie(iedg);
            idneigh = (neigh_ie_unq == neighbor.neigh{ie}(iedg));
            ii_index_neigh{ie_sh}(1:femregion.nbases,:) = repmat(index, 1 ,femregion.nbases);
            jj_index_neigh{ie_sh}(1:femregion.nbases,:) = repmat(index',femregion.nbases,1);

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

                Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
                Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

                Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (nx * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqx + ny * par.mu_ave .* gradedgeqy ))' * phiedgeq;
                Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (ny * par.lambda_ave .* gradedgeqx + nx * par.mu_ave .* gradedgeqy ))' * phiedgeq;
                Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (ny * par.mu_ave .* gradedgeqx + nx * par.lambda_ave .* gradedgeqy ))' * phiedgeq;
                Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (nx * par.mu_ave .* gradedgeqx + ny * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqy ))' * phiedgeq;
                 
            % Absorbing boundary faces
            elseif neigh_ie(iedg) == -3 

                Edge.ABC_S1_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_S1_loc{ie_sh}(1:femregion.nbases,:) + (ds .* par.c11 .* phiedgeq)' * phiedgeq;
                Edge.ABC_S2_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_S2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* par.c12 .* phiedgeq)' * phiedgeq;
                Edge.ABC_S3_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_S3_loc{ie_sh}(1:femregion.nbases,:) + (ds .* par.c21 .* phiedgeq)' * phiedgeq;
                Edge.ABC_S4_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_S4_loc{ie_sh}(1:femregion.nbases,:) + (ds .* par.c22 .* phiedgeq)' * phiedgeq;

                Edge.ABC_R1_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_R1_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.c3_11 .* gradedgeqx + par.c4_11 .* gradedgeqy))' * phiedgeq;
                Edge.ABC_R2_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_R2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.c3_12 .* gradedgeqx + par.c4_12 .* gradedgeqy))' * phiedgeq;
                Edge.ABC_R3_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_R3_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.c3_21 .* gradedgeqx + par.c4_21 .* gradedgeqy))' * phiedgeq;
                Edge.ABC_R4_loc{ie_sh}(1:femregion.nbases,:) = Edge.ABC_R4_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.c3_22 .* gradedgeqx + par.c4_22 .* gradedgeqy))' * phiedgeq;
                                             
            % Elastic neighbor
            elseif neigh_ie(iedg)>0 && femregion.tag(neigh_ie(iedg)) == 'E'
                
                % Element itself
                Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
                Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

                Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (nx * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqx + ny * par.mu_ave .* gradedgeqy ))' * phiedgeq;
                Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (ny * par.lambda_ave .* gradedgeqx + nx * par.mu_ave .* gradedgeqy ))' * phiedgeq;
                Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (ny * par.mu_ave .* gradedgeqx + nx * par.lambda_ave .* gradedgeqy ))' * phiedgeq;
                Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (nx * par.mu_ave .* gradedgeqx + ny * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqy ))' * phiedgeq;

                % Neighboring element
                neigh_idx = find(idneigh)*femregion.nbases+1:(find(idneigh)+1)*(femregion.nbases);
                index_neigh = (neighbor.neigh{ie}(iedg)-1-nel_sh)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
                ii_index_neigh{ie_sh}(neigh_idx,:) = repmat(index, 1,femregion.nbases);
                jj_index_neigh{ie_sh}(neigh_idx,:) = repmat(index_neigh',femregion.nbases,1);

                % Construction of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);
                
                % Neighboring element
                Edge.S1_P_loc{ie_sh}(neigh_idx,:) = Edge.S1_P_loc{ie_sh}(neigh_idx,:) - (ds .* (par.harm_ave .* penalty_geom(iedg)) .* phiedgeq)' * phiedgeqneigh;
                Edge.S4_P_loc{ie_sh}(neigh_idx,:) = Edge.S4_P_loc{ie_sh}(neigh_idx,:) - (ds .* (par.harm_ave .* penalty_geom(iedg)) .* phiedgeq)' * phiedgeqneigh;

                Edge.IT1_loc{ie_sh}(neigh_idx,:) = Edge.IT1_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (nx * (par.lambda_ave + 2*par.mu_ave).* gradedgeqx + ny * par.mu_ave .* gradedgeqy ))' * phiedgeqneigh;
                Edge.IT2_loc{ie_sh}(neigh_idx,:) = Edge.IT2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (ny * par.lambda_ave .* gradedgeqx + nx * par.mu_ave .* gradedgeqy ))' * phiedgeqneigh;
                Edge.IT3_loc{ie_sh}(neigh_idx,:) = Edge.IT3_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (ny * par.mu_ave .* gradedgeqx + nx * par.lambda_ave .* gradedgeqy ))' * phiedgeqneigh;
                Edge.IT4_loc{ie_sh}(neigh_idx,:) = Edge.IT4_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (nx * par.mu_ave .* gradedgeqx + ny * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqy ))' * phiedgeqneigh;
              
            elseif neigh_ie(iedg)>0 && femregion.tag(neigh_ie(iedg)) == 'P'
                
                % Element itself
                Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
                Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

                Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (nx * (par.lambda_ave + 2*par.mu_ave) .* phiedgeq))' *gradedgeqx ...
                                                                                                      + (ds .* (ny * par.mu_ave) .* phiedgeq)' * gradedgeqy;
                Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (ny * par.lambda_ave .* phiedgeq))' * gradedgeqx ...
                                                                                                      + (ds .* (nx * par.mu_ave .* phiedgeq))' *gradedgeqy;
                Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (ny * par.mu_ave .* phiedgeq))' *gradedgeqx ...
                                                                                                      + (ds .* (nx * par.lambda_ave .* phiedgeq))' * gradedgeqy;
                Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (nx * par.mu_ave .*  phiedgeq))' * gradedgeqx ...
                                                                                                      + (ds .* (ny * (par.lambda_ave + 2*par.mu_ave) .* phiedgeq))' * gradedgeqy;

            end
        end

    end
    
end


ii_index  = reshape(cell2mat(ii_index),[femregion.nbases,femregion.nbases*femregion.nel_e]);
jj_index  = reshape(cell2mat(jj_index),[femregion.nbases,femregion.nbases*femregion.nel_e]);

El.V1_loc = reshape(cell2mat(El.V1_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
El.V2_loc = reshape(cell2mat(El.V2_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
El.V3_loc = reshape(cell2mat(El.V3_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
El.V4_loc = reshape(cell2mat(El.V4_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);

El.M1_P_rho_loc  = reshape(cell2mat(El.M1_P_rho_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
El.MPrjP_1_loc   = reshape(cell2mat(El.MPrjP_1_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
                  
El.D1_loc  = reshape(cell2mat(El.D1_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
El.C1_loc  = reshape(cell2mat(El.C1_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);

Edge.ABC_R1_loc = reshape(cell2mat(Edge.ABC_R1_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
Edge.ABC_R2_loc = reshape(cell2mat(Edge.ABC_R2_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
Edge.ABC_R3_loc = reshape(cell2mat(Edge.ABC_R3_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
Edge.ABC_R4_loc = reshape(cell2mat(Edge.ABC_R4_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);

Edge.ABC_S1_loc = reshape(cell2mat(Edge.ABC_S1_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
Edge.ABC_S2_loc = reshape(cell2mat(Edge.ABC_S2_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
Edge.ABC_S3_loc = reshape(cell2mat(Edge.ABC_S3_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);
Edge.ABC_S4_loc = reshape(cell2mat(Edge.ABC_S4_loc),[femregion.nbases,femregion.nbases*femregion.nel_e]);

ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);
jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);

Edge.S1_P_loc  = reshape(cell2mat(Edge.S1_P_loc),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);
Edge.S4_P_loc  = reshape(cell2mat(Edge.S4_P_loc),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);

Edge.IT1_loc  = reshape(cell2mat(Edge.IT1_loc),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);
Edge.IT2_loc  = reshape(cell2mat(Edge.IT2_loc),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);
Edge.IT3_loc  = reshape(cell2mat(Edge.IT3_loc),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);
Edge.IT4_loc  = reshape(cell2mat(Edge.IT4_loc),[femregion.nbases,femregion.nel_e*max_nedges*femregion.nbases]);

del = all(jj_index_neigh == 0,1);
ii_index_neigh(:,del) = [];
jj_index_neigh(:,del) = [];

Edge.S1_P_loc(:,del) = [];
Edge.S4_P_loc(:,del) = [];

Edge.IT1_loc(:,del) = [];
Edge.IT2_loc(:,del) = [];
Edge.IT3_loc(:,del) = [];
Edge.IT4_loc(:,del) = [];

A.V1 = sparse(ii_index,jj_index,El.V1_loc,femregion.ndof_e,femregion.ndof_e);
A.V2 = sparse(ii_index,jj_index,El.V2_loc,femregion.ndof_e,femregion.ndof_e);
A.V3 = sparse(ii_index,jj_index,El.V3_loc,femregion.ndof_e,femregion.ndof_e);
A.V4 = sparse(ii_index,jj_index,El.V4_loc,femregion.ndof_e,femregion.ndof_e);

A.M1_P_rho  = sparse(ii_index,jj_index,El.M1_P_rho_loc,femregion.ndof_e,femregion.ndof_e);
A.MPrjP_1   = sparse(ii_index,jj_index,El.MPrjP_1_loc,femregion.ndof_e,femregion.ndof_e);
                  
A.D1  = sparse(ii_index,jj_index,El.D1_loc,femregion.ndof_e,femregion.ndof_e);
A.C1  = sparse(ii_index,jj_index,El.C1_loc,femregion.ndof_e,femregion.ndof_e);

A.ABC_R1   = sparse(ii_index,jj_index,Edge.ABC_R1_loc,femregion.ndof_e,femregion.ndof_e);
A.ABC_R2   = sparse(ii_index,jj_index,Edge.ABC_R2_loc,femregion.ndof_e,femregion.ndof_e);
A.ABC_R3   = sparse(ii_index,jj_index,Edge.ABC_R3_loc,femregion.ndof_e,femregion.ndof_e);
A.ABC_R4   = sparse(ii_index,jj_index,Edge.ABC_R4_loc,femregion.ndof_e,femregion.ndof_e);

A.ABC_S1   = sparse(ii_index,jj_index,Edge.ABC_S1_loc,femregion.ndof_e,femregion.ndof_e);
A.ABC_S2   = sparse(ii_index,jj_index,Edge.ABC_S2_loc,femregion.ndof_e,femregion.ndof_e);
A.ABC_S3   = sparse(ii_index,jj_index,Edge.ABC_S3_loc,femregion.ndof_e,femregion.ndof_e);
A.ABC_S4   = sparse(ii_index,jj_index,Edge.ABC_S4_loc,femregion.ndof_e,femregion.ndof_e);

A.S1_P     = sparse(ii_index_neigh,jj_index_neigh,Edge.S1_P_loc,femregion.ndof_e,femregion.ndof_e);
A.S2_P     = sparse(femregion.ndof_e,femregion.ndof_e);
A.S3_P     = sparse(femregion.ndof_e,femregion.ndof_e);
A.S4_P     = sparse(ii_index_neigh,jj_index_neigh,Edge.S4_P_loc,femregion.ndof_e,femregion.ndof_e);

A.IT1_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT1_loc,femregion.ndof_e,femregion.ndof_e);
A.IT2_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT2_loc,femregion.ndof_e,femregion.ndof_e);
A.IT3_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT3_loc,femregion.ndof_e,femregion.ndof_e);
A.IT4_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT4_loc,femregion.ndof_e,femregion.ndof_e);

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


Matrices = struct( ...
    'A_E', V + S_P - IT_P - transpose(IT_P), ...
    'Dvel', D, ...
    'Ddis', C, ...
    'Svel', ABC_S, ...
    'Rdis', ABC_R, ...
    'MPrjP',MprjP,...
    'DGe', V + S_P,...
    'M_P_rho', M_P_rho);

fprintf('\n');

