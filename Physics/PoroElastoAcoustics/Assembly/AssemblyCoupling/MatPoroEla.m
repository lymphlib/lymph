%> @file  MatPoroEla.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 9 August 2024
%> @brief Assembly of the matrices for the poroelastic-elastic problem
%>
%==========================================================================
%> @section classMatPoroEla Class description
%==========================================================================
%> @brief Assembly of the matrices for the poroelastic-elastic problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Coupling matrices
%>
%==========================================================================

function [Matrices] = MatPoroEla(Data, neighbor, femregion)

%% Quadrature values
[ref_qNodes_1D, w_1D, ~, ~] = Quadrature(femregion.nqn);

    %% number of interface - faces :
    n_interf=0;
    for ie = 1 : femregion.nel
        if femregion.tag(ie) == 'P'
            neigh_ie = neighbor.neigh{ie};
            n_interf = n_interf + sum(neigh_ie(neigh_ie>0)'.*(femregion.tag(neigh_ie(neigh_ie>0)) == 'E')>0);
        end
    end

    %% Initialization of the matrices
    ii_index_neigh = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    jj_index_neigh = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));

    %% Initialization of the matrices
    C1_P_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C2_P_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C3_P_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C4_P_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
        
    C1_E_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C2_E_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C3_E_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C4_E_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));

    %% Loop over the elements

    ie_gamma = 1;
    id_shift = max([Data.TagElAcu,Data.TagElPoro]);
    if(isempty(id_shift)); id_shift = 0; end
    nel_sh = femregion.nel_p + femregion.nel_a;

    % Visualization of computational progress
    prog = 0;
    fprintf(1,'\t Computation Progress: %3d%%\n',prog);
    
    for ie = 1 : femregion.nel

        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);

        % Id and tag selection for elements
        id_ie  = femregion.id(ie);
        tag_ie = femregion.tag(ie);

        % Check if the element is poroelastic
        if tag_ie == 'P'

            % Selection of the matrix positions associated to element ie
            index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
            
            % Extraction of neighbor element and their edges
            neigh_ie      = neighbor.neigh{ie};
            neighedges_ie = neighbor.neighedges{ie};

            % Extraction of element geometrical information
            coords_ie          = femregion.coords_element{ie};
            [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);
   
            %% Coupling integrals
            
            % Computation of all the penalty coefficients for the element ie
            [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

            for iedg = 1 : neighbor.nedges(ie) % loop over faces

                % Extraction of tag and id of neighbor el
                id_el_neigh = neigh_ie(iedg);

                if id_el_neigh > 0
                    id_neigh  = femregion.id(id_el_neigh);
                    tag_neigh = femregion.tag(id_el_neigh);
                else
                    id_neigh  = 0;
                    tag_neigh = 'NaN';
                end
    
                if neigh_ie(iedg)>0 && tag_neigh == 'E'
    
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
                    par.mu     = Data.mu{id_ie}(xq,yq);
                    par.lam    = Data.lam{id_ie}(xq,yq);
                    par.mu_n   = Data.mu_el{id_neigh-id_shift}(xq,yq);
                    par.lam_n  = Data.lam_el{id_neigh-id_shift}(xq,yq);
    
                    % Auxiliary quantities (cf. physical parameters) for vector assembling
                    par.lam_ave  = 2 * par.lam .* par.lam_n ./ (par.lam + par.lam_n);
                    par.mu_ave   = 2 * par.mu .* par.mu_n ./ (par.mu + par.mu_n);
                    par.harm_ave = par.lam_ave + 2*par.mu_ave;
    
                    % Construction and evalutation on the quadrature points of the basis functions
                    phiedgeq = Evalshape2D(femregion, ie, qNodes_1D);

                    % Neighboring element
                    index_neigh = (neighbor.neigh{ie}(iedg)-1-nel_sh)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

                    % Extraction of the indexes for assembling face matrices
                    ii_index_neigh{ie_gamma} = repmat(index, 1,femregion.nbases);
                    jj_index_neigh{ie_gamma} = repmat(index_neigh',femregion.nbases,1);

                    % Construction and evalutation on the quadrature points of the basis functions for the neighbor
                    [phiedgeqneigh, gradedgeneighqx, gradedgeneighqy] = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);
            
                    C1_P_loc{ie_gamma} = - ( ds .* (par.lam_ave + 2*par.mu_ave) .* nx .* phiedgeq)' * gradedgeneighqx ...
                                         - ( ds .* par.mu_ave .* ny .* phiedgeq)' * gradedgeneighqy ...
                                         - penalty_geom(iedg) * ( ds .* par.harm_ave .* phiedgeq)' * phiedgeqneigh;
    
                    C2_P_loc{ie_gamma} = - ( ds .* par.lam_ave .* nx .* phiedgeq)' * gradedgeneighqy ...
                                         - ( ds .* par.mu_ave .* ny .* phiedgeq)' * gradedgeneighqx;
    
                    C3_P_loc{ie_gamma} = - ( ds .* par.mu_ave .* nx .* phiedgeq)' * gradedgeneighqy ...
                                         - ( ds .* par.lam_ave .* ny .* phiedgeq)' * gradedgeneighqx;
    
                    C4_P_loc{ie_gamma} = - ( ds .* par.mu_ave .* nx .* phiedgeq)' * gradedgeneighqx ...
                                         - ( ds .* (par.lam_ave + 2*par.mu_ave) .* ny .* phiedgeq)' * gradedgeneighqy ...
                                         - penalty_geom(iedg) * ( ds .* par.harm_ave .* phiedgeq)' * phiedgeqneigh;

                    C1_E_loc{ie_gamma} = - ( ds .* (par.lam_ave + 2*par.mu_ave) .* nx .* gradedgeneighqx)' * phiedgeq ...
                                         - ( ds .* par.mu_ave .* ny .* gradedgeneighqy)' *  phiedgeq ...
                                         - penalty_geom(iedg) * ( ds .* par.harm_ave .* phiedgeqneigh)' * phiedgeq;
    
                    C2_E_loc{ie_gamma} = - ( ds .* par.mu_ave .* nx .* gradedgeneighqy)' * phiedgeq ...
                                         - ( ds .* par.lam_ave .* ny .* gradedgeneighqx)' * phiedgeq;

                    C3_E_loc{ie_gamma} = - ( ds .* par.lam_ave .* nx .* gradedgeneighqy)' * phiedgeq ...
                                         - ( ds .* par.mu_ave .* ny .* gradedgeneighqx)' * phiedgeq;

                    C4_E_loc{ie_gamma} = - ( ds .* par.mu_ave .* nx .* gradedgeneighqx)' * phiedgeq ...
                                         - ( ds .* (par.lam_ave + 2*par.mu_ave) .* ny .* gradedgeneighqy)' *  phiedgeq ...
                                         - penalty_geom(iedg) * ( ds .* par.harm_ave .* phiedgeqneigh)' * phiedgeq;
                    
                    ie_gamma = ie_gamma + 1;
    
                end
            end
        end
    end
   
    % Local matrix to global matrix
    ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,n_interf*femregion.nbases]);
    jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,n_interf*femregion.nbases]);

    C1_P_loc  = reshape(cell2mat(C1_P_loc),[femregion.nbases,n_interf*femregion.nbases]);
    C2_P_loc  = reshape(cell2mat(C2_P_loc),[femregion.nbases,n_interf*femregion.nbases]);
    C3_P_loc  = reshape(cell2mat(C3_P_loc),[femregion.nbases,n_interf*femregion.nbases]);
    C4_P_loc  = reshape(cell2mat(C4_P_loc),[femregion.nbases,n_interf*femregion.nbases]);

    C1_E_loc  = reshape(pagetranspose(cell2mat(C1_E_loc)),[femregion.nbases,n_interf*femregion.nbases]);
    C2_E_loc  = reshape(pagetranspose(cell2mat(C2_E_loc)),[femregion.nbases,n_interf*femregion.nbases]);
    C3_E_loc  = reshape(pagetranspose(cell2mat(C3_E_loc)),[femregion.nbases,n_interf*femregion.nbases]);
    C4_E_loc  = reshape(pagetranspose(cell2mat(C4_E_loc)),[femregion.nbases,n_interf*femregion.nbases]);

    C1_P = sparse(ii_index_neigh,jj_index_neigh,C1_P_loc,femregion.ndof_p,femregion.ndof_e);
    C2_P = sparse(ii_index_neigh,jj_index_neigh,C2_P_loc,femregion.ndof_p,femregion.ndof_e);
    C3_P = sparse(ii_index_neigh,jj_index_neigh,C3_P_loc,femregion.ndof_p,femregion.ndof_e);
    C4_P = sparse(ii_index_neigh,jj_index_neigh,C4_P_loc,femregion.ndof_p,femregion.ndof_e);

    C1_E = sparse(jj_index_neigh,ii_index_neigh,C1_E_loc,femregion.ndof_e,femregion.ndof_p);
    C2_E = sparse(jj_index_neigh,ii_index_neigh,C2_E_loc,femregion.ndof_e,femregion.ndof_p);
    C3_E = sparse(jj_index_neigh,ii_index_neigh,C3_E_loc,femregion.ndof_e,femregion.ndof_p);
    C4_E = sparse(jj_index_neigh,ii_index_neigh,C4_E_loc,femregion.ndof_e,femregion.ndof_p);
    
    %% Building the matrices
    
    C_P = [C1_P  C2_P; C3_P C4_P];
    C_E = [C1_E  C2_E; C3_E C4_E];
    
    Matrices = struct('C_P', C_P,...
        'C_E', C_E);
    
    fprintf('\n');
    
    
    
    
    
    

