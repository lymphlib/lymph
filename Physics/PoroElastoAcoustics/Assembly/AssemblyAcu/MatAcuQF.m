%> @file  MatAcuQF.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 13 August 2024
%> @brief Assembly of the matrices for the acoustic problem \cite ABNM2021
%>
%==========================================================================
%> @section classMatAcu Class description
%==========================================================================
%> @brief Assembly of the matrices for the acoustic problem \cite ABNM2021
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Struct Matrices computed for the poro domain
%
%> @retval Matrices  Acoustic matrices (Mass, Stiffness, dG, projection matrix, etc)
%>
%==========================================================================

function [Matrices] = MatAcuQF(Data, neighbor, femregion)

%% 1D Quadrature values
[ref_qNodes_1D, w_1D, ~, ~] = Quadrature(femregion.nqn);

%% Construction of the basis functions
[Lx, Ly, dLx, dLy] = Evalshape2DCoeff(femregion.degree);

%% Construction of coefficients of matrices terms
[Coeff] = MatricesCoeff(femregion, Lx, Ly, dLx, dLy);

%% Initialization of the matrices
max_nedges = max(neighbor.nedges+1);

ii_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));
jj_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_a),femregion.nbases,femregion.nbases,ones(1,femregion.nel_a));
ii_index_QF = zeros(femregion.nbases^2,femregion.nel_a);
jj_index_QF = zeros(femregion.nbases^2,femregion.nel_a);
ii_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_a),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_a));
jj_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_a),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_a));

%% Initialization of the matrices
El.M_A_loc   = zeros(femregion.nbases^2,femregion.nel_a);
El.MPrjA_loc = zeros(femregion.nbases^2,femregion.nel_a);

El.W_loc     = zeros(femregion.nbases^2,femregion.nel_a);

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
        ii_index_QF(:,ie_sh) = reshape(repmat(index, 1,femregion.nbases),[femregion.nbases^2 1]);
        jj_index_QF(:,ie_sh) = reshape(repmat(index',femregion.nbases,1),[femregion.nbases^2 1]);
    

        % Extraction of neighbor element and their edges
        neigh_ie      = neighbor.neigh{ie};
        neigh_ie_unq  = unique(neighbor.neigh{ie});
        neighedges_ie = neighbor.neighedges{ie};

        % Construction of the polygon in reference coordinates
        [coords_ie, Jdet, BJ, ~] = GetJacobianPolygon(femregion.bbox(ie,:), femregion.coords_element{ie}');

        % Computation of the monomial integrals in the polygon v
        [I, ~] = IntegralOverPolygon(1, 2*femregion.degree, 2*femregion.degree, coords_ie);

        % Computation of the bases integrals
        Integral.phiphiC = Jdet*Coeff.phiphiC*I;
        Integral.gradxgradxC = Jdet/(BJ(1,1))^2*Coeff.gradxgradxC*I;
        Integral.gradxgradyC = Jdet/(BJ(1,1)*BJ(2,2))*Coeff.gradxgradyC*I;
        Integral.gradygradxC = Jdet/(BJ(1,1)*BJ(2,2))*Coeff.gradygradxC*I;
        Integral.gradygradyC = Jdet/(BJ(2,2))^2*Coeff.gradygradyC*I;
        Integral.gradxphiC   = Jdet/(BJ(1,1))*Coeff.gradxphiC*I;
        Integral.gradyphiC   = Jdet/(BJ(2,2))*Coeff.gradyphiC*I;

        % Evaluation of physical parameters
        par.rho_a   = Data.rho_a{id_ie-id_shift}(0,0);
        par.c       = Data.c{id_ie-id_shift}(0,0);
                    
        %% Local matrix assembling
        El.W_loc(:,ie_sh) = par.rho_a*(Integral.gradxgradxC+Integral.gradygradyC);

        El.M_A_loc(:,ie_sh)   = (par.c.^(-2).* par.rho_a)*Integral.phiphiC;
        El.MPrjA_loc(:,ie_sh) = Integral.phiphiC;
 
        %% Boundary integrals and stabilization terms
        
        % Extraction of element geometrical information
        coords_ie          = femregion.coords_element{ie};
        [normals,meshsize] = GetNormalsMeshSizeFaces(femregion.coords_element{ie});

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
            par.rho_a   = Data.rho_a{id_ie-id_shift}(xq,yq);
            par.c       = Data.c{id_ie-id_shift}(xq,yq);

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
                        
            % Construction of the basis functions
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
                ii_index_neigh{ie_sh}(neigh_idx,:) = repmat(index, 1,femregion.nbases);
                jj_index_neigh{ie_sh}(neigh_idx,:) = repmat(index_neigh',femregion.nbases,1);

                % Construction of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);

                % Neighboring element
                Edge.S_A_loc{ie_sh}(neigh_idx,:)  = Edge.S_A_loc{ie_sh}(neigh_idx,:) - (ds .* (par.rho_a_ave .* penalty_geom(iedg)) .* phiedgeq)' * phiedgeqneigh; 
                Edge.IT_A_loc{ie_sh}(neigh_idx,:) = Edge.IT_A_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (nx * par.rho_a_ave .* gradedgeqx + ny * par.rho_a_ave .* gradedgeqy))' * phiedgeqneigh;

            end
        end

    end
    
end
                 
ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);
jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);

Edge.S_A_loc  = reshape(cell2mat(Edge.S_A_loc),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);
Edge.IT_A_loc = reshape(cell2mat(Edge.IT_A_loc),[femregion.nbases,femregion.nel_a*max_nedges*femregion.nbases]);

del = all(jj_index_neigh == 0,1);
ii_index_neigh(:,del) = [];
jj_index_neigh(:,del) = [];

Edge.S_A_loc(:,del) = [];
Edge.IT_A_loc(:,del) = [];

A.W = sparse(ii_index_QF,jj_index_QF,El.W_loc,femregion.ndof_a,femregion.ndof_a);

A.M_A   = sparse(ii_index_QF,jj_index_QF,El.M_A_loc,femregion.ndof_a,femregion.ndof_a);
A.MPrjA = sparse(ii_index_QF,jj_index_QF,El.MPrjA_loc,femregion.ndof_a,femregion.ndof_a);
                  
A.S_A   = sparse(ii_index_neigh,jj_index_neigh,Edge.S_A_loc,femregion.ndof_a,femregion.ndof_a);
A.IT_A  = sparse(ii_index_neigh,jj_index_neigh,Edge.IT_A_loc,femregion.ndof_a,femregion.ndof_a);

Matrices = struct( 'MPrjA', A.MPrjA,...
                   'A_A', A.W + A.S_A - A.IT_A - transpose(A.IT_A),...
                   'M_A', A.M_A, ...
                   'DGa', A.W + A.S_A);

fprintf('\n');

