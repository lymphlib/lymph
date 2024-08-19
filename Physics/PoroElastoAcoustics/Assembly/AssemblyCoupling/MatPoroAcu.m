%> @file  MatPoroAcu.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 9 August 2024
%> @brief Assembly of the matrices for the poroelastic-acoustic problem \cite ABNM2021
%>
%==========================================================================
%> @section classMatPoroAcu Class description
%==========================================================================
%> @brief Assembly of the matrices for the poroelastic-acoustic problem \cite ABNM2021
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Struct containing the matrices of the problem
%
%> @retval Matrices  Coupling matrices
%>
%==========================================================================

function [Matrices] = MatPoroAcu(Data, neighbor, femregion)

    %% Quadrature values
    [ref_qNodes_1D, w_1D, ~, ~] = Quadrature(femregion.nqn);

    %% number of interface - faces :
    n_interf=0;
    for ie = 1 : femregion.nel
        if femregion.tag(ie) == 'P'
            neigh_ie = neighbor.neigh{ie};
            n_interf = n_interf + sum(neigh_ie(neigh_ie>0)'.*(femregion.tag(neigh_ie(neigh_ie>0)) == 'A')>0);
        end
    end

    %% Initialization of the matrices
    ii_index_neigh = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    jj_index_neigh = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));

    %% Initialization of the matrices
    C1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
    C2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,n_interf),femregion.nbases,femregion.nbases,ones(1,n_interf));
   
    %% Loop over the elements

    ie_gamma = 1;
    id_shift = max(Data.TagElPoro);
    if(isempty(id_shift)); id_shift = 0; end
    nel_sh = femregion.nel_p;

    % Visualization of computational progress
    prog = 0;
    fprintf(1,'\t Computation Progress: %3d%%\n',prog);

    for ie = 1 : femregion.nel

        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);

        % Id and tag selection for elements
        tag_ie = femregion.tag(ie);

        % Check if the element is poroelastic
        if tag_ie == 'P'

            % Selection of the matrix positions associated to element ie
            index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
            
            % Extraction of neighbor element and their edges
            neigh_ie      = neighbor.neigh{ie};

            % Extraction of element geometrical information
            coords_ie          = femregion.coords_element{ie};
            [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);

            %% Coupling integrals

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

                if neigh_ie(iedg) >0 && tag_neigh == 'A'

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
                    par.rho_a = Data.rho_a{id_neigh-id_shift}(xq,yq);

                    % Construction of the basis functions
                    phiedgeq = Evalshape2D(femregion, ie, qNodes_1D);

                    % Neighboring element
                    index_neigh = (neighbor.neigh{ie}(iedg)-1-nel_sh)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
                    ii_index_neigh{ie_gamma} = repmat(index, 1,femregion.nbases);
                    jj_index_neigh{ie_gamma} = repmat(index_neigh',femregion.nbases,1);

                    % Construction of the basis functions for the neighbor
                    phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);

                    C1_loc{ie_gamma} = (ds .* par.rho_a .* nx .* phiedgeq)' * phiedgeqneigh;
                    C2_loc{ie_gamma} = (ds .* par.rho_a .* ny .* phiedgeq)' * phiedgeqneigh;

                    ie_gamma = ie_gamma + 1;

                end
            end
        end
    end

    %% Reshape of the matrices

    ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,n_interf*femregion.nbases]);
    jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,n_interf*femregion.nbases]);

    C1_loc  = reshape(cell2mat(C1_loc),[femregion.nbases,n_interf*femregion.nbases]);
    C2_loc  = reshape(cell2mat(C2_loc),[femregion.nbases,n_interf*femregion.nbases]);

    C1 = sparse(ii_index_neigh,jj_index_neigh,C1_loc,femregion.ndof_p,femregion.ndof_a);
    C2 = sparse(ii_index_neigh,jj_index_neigh,C2_loc,femregion.ndof_p,femregion.ndof_a);

    %% Construction of coupling matrices

    C1_P = [C1 ; C2];
    C2_P = [C1 ; C2];

    C1_A = - C1_P';
    C2_A = - C2_P';

    Matrices = struct('C1_P', C1_P,...
        'C2_P', C2_P,...
        'C1_A', C1_A,...
        'C2_A', C2_A);

    fprintf('\n');







