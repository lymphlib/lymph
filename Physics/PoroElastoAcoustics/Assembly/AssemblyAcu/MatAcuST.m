%> @file  MatAcuST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 6 August 2024
%> @brief Assembly of the matrices for the acoustic problem \cite ABNM2021
%>
%==========================================================================
%> @section classMatAcuST Class description
%==========================================================================
%> @brief Assembly of the matrices for the acoustic problem \cite ABNM2021
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Acoustic matrices (Mass, Stiffness, dG, projection matrix, etc)
%>
%==========================================================================

function [Matrices] = MatAcuST(Data, neighbor, femregion)

%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices
max_nedges = max(neighbor.nedges+1);

ii_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));
jj_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));
ii_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_a),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_a));
jj_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_a),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_a));

%% Initialization of the matrices
El.M_A_loc   = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));
El.MPrjA_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));

El.W_loc     = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));

Edge.S_A_loc   = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_a),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_a));
Edge.IT_A_loc  = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_a),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_a));
 
%% Loop over the elements

% Visualization of computational progress
index_shift = 0;
id_shift = max([Data.TagElPoro]);
if(isempty(id_shift)); id_shift = 0; end
nel_sh = femregion.nel_p;

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
    if tag_ie == 'A'

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
            par.rho_a   = Data.rho_a{id_ie-id_shift}(xq,yq);
            par.c       = Data.c{id_ie-id_shift}(xq,yq);
            
            % Construction and evalutation on the quadrature points of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
                    
            %% Local matrix assembling
            El.W_loc{ie_sh} = El.W_loc{ie_sh} + (dx .* (par.rho_a .* gradqx))' * gradqx + (dx .* (par.rho_a .* gradqy))' * gradqy;
            
            El.M_A_loc{ie_sh}   = El.M_A_loc{ie_sh} + (dx .* (par.c.^(-2).* par.rho_a  .* phiq))' * phiq;
            El.MPrjA_loc{ie_sh} = El.MPrjA_loc{ie_sh} + (dx .* phiq)' * phiq;
            
        end
 
        %% Boundary integrals and stabilization terms
        
        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

        for iedg = 1 : neighbor.nedges(ie) % loop over faces
            
            % Extraction of tag and id of neighbor el
            id_el_neigh = neigh_ie(iedg);
            idneigh = (neigh_ie_unq == neighbor.neigh{ie}(iedg));

            % Extraction of the indexes for assembling face matrices (contribution of the element ie)
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
            par.rho_a   = Data.rho_a{id_ie-id_shift}(xq,yq);
            par.c       = Data.c{id_ie-id_shift}(xq,yq);

            % Auxiliary quantities (cf. physical parameters) for vector assembling
            if strcmp(tag_neigh,'E') %elastic neighbor
                par.rho_a_n   = Data.rho_a{id_ie-id_shift}(xq,yq);
                par.c_n       = Data.c{id_ie-id_shift}(xq,yq);
            elseif strcmp(tag_neigh,'P') %poroelastic neighbor
                par.rho_a_n   = Data.rho_a{id_ie-id_shift}(xq,yq);
                par.c_n       = Data.c{id_ie-id_shift}(xq,yq);
            elseif strcmp(tag_neigh,'A') %acoustic edge
                par.rho_a_n   = Data.rho_a{id_neigh-id_shift}(xq,yq);
                par.c_n       = Data.c{id_neigh-id_shift}(xq,yq);
            elseif strcmp(tag_neigh,'NaN') %boundary edge
                par.rho_a_n   = Data.rho_a{id_ie-id_shift}(xq,yq);
                par.c_n       = Data.c{id_ie-id_shift}(xq,yq);
            end
            par.rho_a_ave = 2*par.rho_a.*par.rho_a_n./(par.rho_a + par.rho_a_n);
            par.rho_c_ave = 2*par.c.*par.c_n./(par.c + par.c_n);      
                        
            % Construction and evalutation on the quadrature points of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);            
            
            % Dirichlet boundary faces
            if neigh_ie(iedg) == -1

                Edge.S_A_loc{ie_sh}(1:femregion.nbases,:)  = Edge.S_A_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.rho_a_ave .* phiedgeq)' * phiedgeq; 
                Edge.IT_A_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT_A_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (nx * par.rho_a_ave .* gradedgeqx + ny * par.rho_a_ave .* gradedgeqy))' * phiedgeq;
            
            % Absorbing boundary faces
            elseif neigh_ie(iedg) == -3 

                disp('Abosrbing conditions not implemented for the acoustic domain')

            % Elastic neighbor
            elseif neigh_ie(iedg)>0 && femregion.tag(neigh_ie(iedg)) == 'A'
                
                % Element itself
                Edge.S_A_loc{ie_sh}(1:femregion.nbases,:)  = Edge.S_A_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.rho_a_ave .* phiedgeq)' * phiedgeq; 
                Edge.IT_A_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT_A_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (nx * par.rho_a_ave .* gradedgeqx + ny * par.rho_a_ave .* gradedgeqy))' * phiedgeq;
   
                % Neighboring element
                neigh_idx = find(idneigh)*femregion.nbases+1:(find(idneigh)+1)*(femregion.nbases);
                index_neigh = (neighbor.neigh{ie}(iedg)-1-nel_sh)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

                % Extraction of the indexes for assembling face matrices (contribution of the neighboring element)
                ii_index_neigh{ie_sh}(neigh_idx,:) = repmat(index, 1,femregion.nbases);
                jj_index_neigh{ie_sh}(neigh_idx,:) = repmat(index_neigh',femregion.nbases,1);

                % Construction and evalutation on the quadrature points of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);

                % Neighboring element
                Edge.S_A_loc{ie_sh}(neigh_idx,:)  = Edge.S_A_loc{ie_sh}(neigh_idx,:) - (ds .* (par.rho_a_ave .* penalty_geom(iedg)) .* phiedgeq)' * phiedgeqneigh; 
                Edge.IT_A_loc{ie_sh}(neigh_idx,:) = Edge.IT_A_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (nx * par.rho_a_ave .* gradedgeqx + ny * par.rho_a_ave .* gradedgeqy))' * phiedgeqneigh;

            end
        end

    end
    
end

% Local matrix to global matrix
ii_index  = reshape(cell2mat(ii_index),[femregion.nbases,femregion.nbases*femregion.nel_a]);
jj_index  = reshape(cell2mat(jj_index),[femregion.nbases,femregion.nbases*femregion.nel_a]);

El.W_loc = reshape(cell2mat(El.W_loc),[femregion.nbases,femregion.nbases*femregion.nel_a]);

El.M_A_loc  = reshape(cell2mat(El.M_A_loc),[femregion.nbases,femregion.nbases*femregion.nel_a]);
El.MPrjA_loc   = reshape(cell2mat(El.MPrjA_loc),[femregion.nbases,femregion.nbases*femregion.nel_a]);
                  
ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);
jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);

Edge.S_A_loc  = reshape(cell2mat(Edge.S_A_loc),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);
Edge.IT_A_loc = reshape(cell2mat(Edge.IT_A_loc),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);

del = all(jj_index_neigh == 0,1);
ii_index_neigh(:,del) = [];
jj_index_neigh(:,del) = [];

Edge.S_A_loc(:,del) = [];
Edge.IT_A_loc(:,del) = [];

A.W = sparse(ii_index,jj_index,El.W_loc,femregion.ndof_a,femregion.ndof_a);

A.M_A   = sparse(ii_index,jj_index,El.M_A_loc,femregion.ndof_a,femregion.ndof_a);
A.MPrjA = sparse(ii_index,jj_index,El.MPrjA_loc,femregion.ndof_a,femregion.ndof_a);
                  
A.S_A   = sparse(ii_index_neigh,jj_index_neigh,Edge.S_A_loc,femregion.ndof_a,femregion.ndof_a);
A.IT_A  = sparse(ii_index_neigh,jj_index_neigh,Edge.IT_A_loc,femregion.ndof_a,femregion.ndof_a);

Matrices = struct( 'MPrjA', A.MPrjA,...
                   'A_A', A.W + A.S_A - A.IT_A - transpose(A.IT_A),...
                   'M_A', A.M_A, ...
                   'DGa', A.W + A.S_A);

fprintf('\n');

