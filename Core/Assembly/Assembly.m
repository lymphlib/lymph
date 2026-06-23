%> @file  Assembly.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 5 June 2026
%> @brief Assembly of the matrices
%>
%==========================================================================
%> @section classAssembly Class description
%==========================================================================
%> @brief            Assembly of matrices with quadrature-free or 
%> subtriangulation implementation for polygons and Gaussian quadrature 
%> for interface integrals.
%>
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure (see example in MatrixAssemblyLaplacian.m). It must contain:
%> - quadrature: quadrature method to be used for volume assembly (QF/ST)
%> - assemblyvolume: if true the assembly of volume integrals is performed
%> - assemblyfaces: if true the assembly of faces integrals is performed
%> - assemblyinternalfaces: if true the assembly of internal faces integrals is performed
%> - assemblytrilinearforms: if true the assembly of trilinear forms matrices is performed
%> - computegradients: if true the gradients of the bases are computed
%> - computelaplacian: if true the laplacian of the bases are computed
%> - computefacegradients: if true the gradients of the bases on the faces are computed
%> - ass_vol_vec: vector of elements where to assembly the volume integrals
%> - ass_face_vec: vector of elements where to assembly the face integrals
%> - computeoldbases: if true the basis associated with the degrees of previous p-adaptive step are computed
%> - Matrices_adapt_old: matrices of the previous assembly iteration
%> The structure can contain additional information whenever needed. For
%> example it should contain the current time in forcing term assembly
%> @param Funcs      Struct of specific assembly functions, the required
%> functions are the following:
%>  - Preallocation: this function preallocates the structure Matrices for 
%> the local matrices assembly divided into Matrices.Volume and Matrices.Faces
%> for volume and faces terms, respectively (see detailed example in 
%> MatrixPreallocationLaplacian.m)
%>  - VolumeAssemblyQF: this function assembles in the structure Matrices.Volume the local
%> volume matrices with the quadrature-free approach (see detailed example in
%> VolumeMatricesAssemblyLaplacianQF.m)
%>  - VolumeAssemblyST: this function assembles in the structure Matrices.Volume the local
%> volume matrices with the subtriangulation approach (see detailed example in
%> VolumeMatricesAssemblyLaplacianST.m)
%>  - Volume3LAssemblyQF: this function assembles in the structure Matrices.Volume3L the local
%> volume matrices associated with trilinear forms with the quadrature-free approach
%>  - Volume3LAssemblyST: this function assembles in the structure Matrices.Volume3L the local
%> volume matrices associated with trilinear forms with the subtriangulation approach
%>  - FacesAssembly: this function assembles in the structure Matrices.Faces the local
%> faces matrices (see detailed example in FacesMatricesAssemblyLaplacian.m)
%>  - FinalMatrices: this function takes the Matrices structure
%> and construct in Matrices the global matrices associated to the
%> problem and to be used in the main solver (see detailed example in
%> SIPMatricesLaplacian.m)
%>
%> @retval Matrices  Matrices struct
%>
%==========================================================================

function [Matrices] = Assembly(Data, neighbor, femregion, AssembInfo, Funcs)

    %% Computation of maximum number of bases, degree and edges for preallocation
    max_nbases = max(femregion.nbases);
    max_degree = max(femregion.degree);
    max_nedges = max(neighbor.nedges);
    
    if isfield(femregion,"degree_old")
        degs_unq = unique([femregion.degree;femregion.degree_old])';
    else
        degs_unq = unique(femregion.degree)';
    end
        
    %% Quadrature values
    switch AssembInfo.quadrature

        case "QF"          
            
            % Cell structures Preallocation
            ref_qNodes_1D = cell(max_degree,1);
            w_1D          = cell(max_degree,1);
            Coeff_QF      = cell(max_degree,1);
            
            % Loop over polynomial degrees
            for deg = degs_unq

                % Quadrature nodes computation
                [ref_qNodes_1D{deg}, w_1D{deg}, ~, ~] = Quadrature(deg + 1);

                % Quadrature-free coefficients computation
                [Coeff_QF{deg}] = MakeCoefficientsQF(deg, (deg+1)*(deg+2)/2);
            end
            ref_qNodes_2D = [];
            w_2D = [];

        case "ST"
            
            % Cell structures Preallocation
            ref_qNodes_1D = cell(max_degree,1);
            ref_qNodes_2D = cell(max_degree,1);
            w_1D          = cell(max_degree,1);
            w_2D          = cell(max_degree,1);
            
            % Loop over polynomial degrees
            for deg = degs_unq

                % Quadrature nodes computation
                [ref_qNodes_1D{deg}, w_1D{deg}, ref_qNodes_2D{deg}, w_2D{deg}] = Quadrature(deg + 1);
            end
            Coeff_QF = [];
    end

    %% Initialization of the general functions for matrices and forcing terms assembly

    % General matrices struct containing preallocation for volume and faces terms
    switch AssembInfo.quadrature
        case "QF"
            GenMatrices.VolMatrix = zeros(max_nbases^2,1);
        case "ST"
            GenMatrices.VolMatrix = zeros(max_nbases,max_nbases);
    end

    GenMatrices.FaceMatrix   = zeros(max_nbases*(max_nedges+1),max_nbases);
    GenMatrices.FaceMatrixBd = zeros(max_nbases,max_nbases);
    GenMatrices.Vector       = zeros(max_nbases,1);
    GenMatrices.CellVector   = zeros(1,1);
    GenMatrices.VolMatrix3L  = zeros(max_nbases^2,max_nbases);

    [Matrices] = Funcs.Preallocation(GenMatrices);

    Matrices = repmat({Matrices}, femregion.nel, 1);

    % Indices cells preallocation
    ii_index = repmat({GenMatrices.VolMatrix}, femregion.nel, 1);
    jj_index = repmat({GenMatrices.VolMatrix}, femregion.nel, 1);

    ii_indexBd = repmat({GenMatrices.FaceMatrixBd}, femregion.nel, 1);
    jj_indexBd = repmat({GenMatrices.FaceMatrixBd}, femregion.nel, 1);

    ii_indexVec = repmat({GenMatrices.Vector}, femregion.nel, 1);

    ii_index_neigh = repmat({GenMatrices.FaceMatrix}, femregion.nel, 1);
    jj_index_neigh = repmat({GenMatrices.FaceMatrix}, femregion.nel, 1);

    %% Loop over the elements
    parfor ie = 1:femregion.nel

        % Selection of the matrix positions associated to element ie
        elem = struct();
        elem.index = (sum(femregion.nbases(1:ie-1)))* ones(femregion.nbases(ie),1) + (1:femregion.nbases(ie))';

        switch AssembInfo.quadrature
            case "QF"
                ii_index{ie}(1:femregion.nbases(ie)^2,1) = reshape(repmat(elem.index, 1,femregion.nbases(ie)),[femregion.nbases(ie)^2 1]);
                jj_index{ie}(1:femregion.nbases(ie)^2,1) = reshape(repmat(elem.index',femregion.nbases(ie),1),[femregion.nbases(ie)^2 1]);
            case "ST"
                ii_index{ie}(1:femregion.nbases(ie),1:femregion.nbases(ie)) = repmat(elem.index, 1, femregion.nbases(ie));
                jj_index{ie}(1:femregion.nbases(ie),1:femregion.nbases(ie)) = repmat(elem.index', femregion.nbases(ie), 1);
        end

        ii_indexBd{ie}(1:femregion.nbases(ie),1:femregion.nbases(ie)) = repmat(elem.index, 1, femregion.nbases(ie));
        jj_indexBd{ie}(1:femregion.nbases(ie),1:femregion.nbases(ie)) = repmat(elem.index', femregion.nbases(ie), 1);

        ii_indexVec{ie}(1:femregion.nbases(ie),1) = elem.index;

        % Extraction of neighbor element and their edges
        neigh_ie      = neighbor.neigh{ie};
        neigh_ie_unq  = unique(neighbor.neigh{ie}','rows')';
        neighedges_ie = neighbor.neighedges{ie};

        if AssembInfo.assemblyvolume

            % Control if is needed to assemble or an old assembly can be used to save computational time
            if AssembInfo.ass_vol_vec(ie)

                % Assembly of the volume terms
                switch AssembInfo.quadrature
                    case "QF"

                        % Construction of the bases integrals on the polygonal element
                        [Integral] = MakeIntegralsQF(Coeff_QF{femregion.degree(ie)}, femregion, ie);

                        % Computation of the bilinear form matrices related to volume integrals
                        [Matrices{ie}.Volume] = Funcs.VolumeAssemblyQF(Data, Integral, Matrices{ie}.Volume, ie, femregion.id_phys(ie), femregion.nbases(ie), AssembInfo);

                        if AssembInfo.assemblytrilinearforms
                            % Computation of the trilinear form matrices related to volume integrals
                            [Matrices{ie}.Volume3L] = Funcs.Volume3LAssemblyQF(Data, Integral, Matrices{ie}.Volume3L, ie, femregion.id_phys(ie), femregion.nbases(ie), AssembInfo);
                        end

                    case "ST"

                        % Extraction of element geometrical information
                        coords_ie = femregion.coords_element{ie};

                        % Creation of the subtriangulation of the element
                        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
                        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
                        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

                        % Construction of Jacobian and quadrature nodes
                        [detBJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie, Tria, ref_qNodes_2D{femregion.degree(ie)});

                        elem.xq  = qNodes_2D(:,1);
                        elem.yq  = qNodes_2D(:,2);

                        % Scaled weights
                        elem.dx = reshape(w_2D{femregion.degree(ie)}*detBJ,length(w_2D{femregion.degree(ie)})*length(detBJ),1);

                        % Construction of the basis functions
                        if AssembInfo.computegradients

                            if AssembInfo.computelaplacian
                                % Compute bases, gradients and laplacian
                                [elem.phiq, elem.gradqx, elem.gradqy, elem.lapqxx, elem.lapqxy, elem.lapqyx, elem.lapqyy] = Evalshape2D(femregion, ie, qNodes_2D);
                            else
                                % Compute bases and gradients
                                [elem.phiq, elem.gradqx, elem.gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
                            end

                        else

                            % Compute bases
                            [elem.phiq] = Evalshape2D(femregion, ie, qNodes_2D);

                        end

                        if isfield(AssembInfo,"computeoldbases") && AssembInfo.computeoldbases

                            % Construct the old indices map
                            elem.index_old = (sum(femregion.nbases_old(1:ie-1)))* ones(femregion.nbases_old(ie),1) + (1:femregion.nbases_old(ie))';

                            if femregion.degree_old(ie) > femregion.degree(ie)

                                % Construction of Jacobian and quadrature nodes
                                [detBJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie, Tria, ref_qNodes_2D{femregion.degree_old(ie)});

                                elem.xq  = qNodes_2D(:,1);
                                elem.yq  = qNodes_2D(:,2);

                                % Scaled weights
                                elem.dx = reshape(w_2D{femregion.degree_old(ie)}*detBJ,length(w_2D{femregion.degree_old(ie)})*length(detBJ),1);

                            end

                            femregion_app = femregion;
                            femregion_app.degree = femregion_app.degree_old;

                            % Compute old bases for projection matrix
                            [elem.phiq]     = Evalshape2D(femregion, ie, qNodes_2D);
                            [elem.phiq_old] = Evalshape2D(femregion_app, ie, qNodes_2D);
                        end

                        % Computation of the matrices related to volume integrals
                        [Matrices{ie}.Volume] = Funcs.VolumeAssemblyST(Data, Matrices{ie}.Volume, elem, ie, femregion.id_phys(ie), femregion.nbases(ie), AssembInfo);

                        if AssembInfo.assemblytrilinearforms
                            % Computation of the trilinear form matrices related to volume integrals
                            [Matrices{ie}.Volume3L] = Funcs.Volume3LAssemblyST(Data, Matrices{ie}.Volume3L, elem, ie, femregion.id_phys(ie), femregion.nbases(ie), AssembInfo);
                        end

                end

            else

                % If assembly is not needed construct the old indices map
                elem.index_old = (sum(femregion.nbases_old(1:ie-1)))* ones(femregion.nbases_old(ie),1) + (1:femregion.nbases_old(ie))';

                % Cutting using the nested bases
                if femregion.nbases_old(ie) > femregion.nbases(ie)

                    % x-index vector
                    idx1 = reshape(tril((1:femregion.degree(ie)+2).*ones(femregion.degree(ie)+2,femregion.degree(ie)+2)),[1 (femregion.degree(ie)+2)^2]);
                    idx1(idx1==0) = [];

                    % y-index vector
                    idx2 = reshape(flip(triu((femregion.degree(ie)+2:-1:1).*ones(femregion.degree(ie)+2,femregion.degree(ie)+2)),2)',[1 (femregion.degree(ie)+2)^2]);
                    idx2(idx2==0) = [];

                    filt = (idx1+idx2-2<=femregion.degree(ie));

                    elem.index_old = elem.index_old(filt);

                end

                fields = fieldnames(Matrices{ie}.Volume);

                % Reuse the local volume matrices
                switch AssembInfo.quadrature
                    case "QF"
                        for k = 1:numel(fields)
                            Matrices{ie}.Volume.(fields{k})(1:femregion.nbases(ie)^2,1) = reshape(AssembInfo.Matrices_adapt_old.Volume.(fields{k})(elem.index_old,elem.index_old),femregion.nbases(ie)^2,1);
                        end

                    case "ST"
                        for k = 1:numel(fields)
                            if size(Matrices{ie}.Volume.(fields{k}),2) == 1
                                if size(Matrices{ie}.Volume.(fields{k}),1) == 1
                                    Matrices{ie}.Volume.(fields{k}) = AssembInfo.Matrices_adapt_old.Volume.(fields{k})(ie);
                                else
                                    Matrices{ie}.Volume.(fields{k})(1:femregion.nbases(ie),1) = AssembInfo.Matrices_adapt_old.Volume.(fields{k})(elem.index_old,1);
                                end
                            else
                                Matrices{ie}.Volume.(fields{k})(1:femregion.nbases(ie),1:femregion.nbases(ie)) = AssembInfo.Matrices_adapt_old.Volume.(fields{k})(elem.index_old,elem.index_old);
                            end
                        end
                end

                if AssembInfo.assemblytrilinearforms
                    elem.index_old1 = elem.index_old-sum(femregion.nbases_old(1:ie-1));
                    elem.index_old2 = reshape(max(femregion.nbases_old)*(elem.index_old1'-1)+elem.index_old1,1,[]);

                    fields3L = fieldnames(Matrices{ie}.Volume3L);
                    for k = 1:numel(fields3L)
                        Matrices{ie}.Volume3L.(fields3L{k})(1:femregion.nbases(ie)^2,1:femregion.nbases(ie)) = AssembInfo.Matrices_adapt_old.Volume3L.(fields3L{k}){ie}(elem.index_old2,elem.index_old1);
                    end
                end

            end
        end

        %% Boundary integrals and stabilization terms
        if AssembInfo.assemblyfaces

            face = struct();
            face.ie = ie;

            % Extraction of element geometrical information
            coords_ie          = femregion.coords_element{ie};
            [normals,meshsize] = GetNormalsMeshSizeFaces(femregion.coords_element{ie});

            % Computation of all the penalty coefficients for the element ie
            [face.penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie(1,:), meshsize);

            diagPass = 0;

            % Loop over faces
            for iedg = 1 : neighbor.nedges(ie)

                if neigh_ie(1,iedg)<0 || AssembInfo.assemblyinternalfaces

                    face.iedg = iedg;
                    face.neigh_ie = neigh_ie(iedg);

                    % Extraction of the id of the neighboring element in the neighbor matrices
                    face.idneigh = (neigh_ie_unq(1,:) == neighbor.neigh{ie}(1,iedg));
                    ii_index_neigh{ie}(1:femregion.nbases(ie),1:femregion.nbases(ie)) = repmat(elem.index, 1, femregion.nbases(ie));
                    jj_index_neigh{ie}(1:femregion.nbases(ie),1:femregion.nbases(ie)) = repmat(elem.index', femregion.nbases(ie), 1);

                    % Construction of the indices for the neighbor
                    if neigh_ie(1,iedg) > 0

                        first_index = sum(femregion.nbases(1:neigh_ie(iedg)-1));
                        face.neigh_idx = find(face.idneigh)*femregion.nbases(ie)+1:(find(face.idneigh)+1)*(femregion.nbases(ie));
                        index_neigh = first_index*ones(femregion.nbases(neigh_ie(iedg)),1) + (1:femregion.nbases(neigh_ie(iedg)))';

                        ii_index_neigh{ie}(face.neigh_idx,1:femregion.nbases(neigh_ie(iedg))) = repmat(elem.index, 1, femregion.nbases(neigh_ie(iedg)));
                        jj_index_neigh{ie}(face.neigh_idx,1:femregion.nbases(neigh_ie(iedg))) = repmat(index_neigh', femregion.nbases(ie), 1);
                    end

                    if AssembInfo.ass_face_vec(ie)

                        % Extraction of the edge coordinates
                        if iedg == neighbor.nedges(ie)
                            p1 = coords_ie(iedg,:);
                            p2 = coords_ie(1,:);
                        else
                            p1 = coords_ie(iedg,:);
                            p2 = coords_ie(iedg+1,:);
                        end

                        % Construction of quadrature nodes on the face with the maximum degree between the element and the neighbor
                        if neigh_ie(iedg) > 0 
                            degree_q1D = max(femregion.degree(neigh_ie(iedg)),femregion.degree(ie));
                        else
                            degree_q1D = femregion.degree(ie);
                        end

                        [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D{degree_q1D});

                        face.xq = qNodes_1D(:,1);
                        face.yq = qNodes_1D(:,2);

                        % Scaled weights
                        face.ds = meshsize(iedg) * w_1D{degree_q1D};

                        % Construction of the basis functions
                        if AssembInfo.computefacegradients
                            [face.phiedgeq, face.gradedgeqx, face.gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);
                        else
                            [face.phiedgeq] = Evalshape2D(femregion, ie, qNodes_1D);
                        end

                        % Extraction of normals to the face
                        face.nx = normals(1,iedg);
                        face.ny = normals(2,iedg);

                        % Construction of the basis functions for the neighbor
                        if neigh_ie(1,iedg) > 0 

                            if AssembInfo.computefacegradients
                                [face.phiedgeqneigh, face.gradedgeqxneigh, face.gradedgeqyneigh] = Evalshape2D(femregion, neigh_ie(1,iedg), qNodes_1D);
                            else
                                [face.phiedgeqneigh] = Evalshape2D(femregion, neigh_ie(1,iedg), qNodes_1D);
                            end

                        end

                        % Computation of the matrices related to face integrals
                        [Matrices{ie}.Faces] = Funcs.FacesAssembly(Data, femregion, Matrices{ie}.Faces, face, AssembInfo);

                    else

                        index_old_ie = (sum(femregion.nbases_old(1:ie-1)))* ones(femregion.nbases_old(ie),1) + (1:femregion.nbases_old(ie))';

                        fields = fieldnames(Matrices{ie}.Faces);

                        if diagPass == 0
                            for k = 1:numel(fields)
                                if size(Matrices{ie}.Faces.(fields{k}),2) == 1
                                    if size(Matrices{ie}.Faces.(fields{k}),1) == 1
                                        Matrices{ie}.Faces.(fields{k}) = AssembInfo.Matrices_adapt_old.Faces.(fields{k})(ie);
                                    else
                                        Matrices{ie}.Faces.(fields{k})(1:femregion.nbases(ie),1) = AssembInfo.Matrices_adapt_old.Faces.(fields{k})(index_old_ie,1);
                                    end
                                else
                                    Matrices{ie}.Faces.(fields{k})(1:femregion.nbases(ie),1:femregion.nbases(ie)) = AssembInfo.Matrices_adapt_old.Faces.(fields{k})(index_old_ie,index_old_ie);
                                end
                            end
                            diagPass = 1;
                        end

                        if neigh_ie(1,iedg) > 0 
                            index_old_neigh_ie = (sum(femregion.nbases_old(1:neigh_ie(iedg)-1)))*ones(femregion.nbases_old(neigh_ie(iedg)),1) + (1:femregion.nbases_old(neigh_ie(iedg)))';

                            for k = 1:numel(fields)
                                if not(size(Matrices{ie}.Faces.(fields{k}),2) == 1 && size(Matrices{ie}.Faces.(fields{k}),1) == 1)
                                    Matrices{ie}.Faces.(fields{k})(face.neigh_idx, 1:femregion.nbases(face.neigh_ie)) = AssembInfo.Matrices_adapt_old.Faces.(fields{k})(index_old_ie, index_old_neigh_ie);
                                end
                            end
                        end

                    end
                end
            end
        end
    end

    %% Construction of global matrices
    if isfield(AssembInfo,"Matrices_adapt_old")
        AssembInfo = rmfield(AssembInfo,"Matrices_adapt_old");
    end

    % Construction of final indices structure
    ii_index_neigh = reshape(vertcat(ii_index_neigh{:}), [1, femregion.nel*(max_nedges+1)*max_nbases^2]);
    jj_index_neigh = reshape(vertcat(jj_index_neigh{:}), [1, femregion.nel*(max_nedges+1)*max_nbases^2]);

    ii_index = reshape(vertcat(ii_index{:}),[1, max_nbases^2*femregion.nel]);
    jj_index = reshape(vertcat(jj_index{:}),[1, max_nbases^2*femregion.nel]);

    ii_indexBd = reshape(vertcat(ii_indexBd{:}),[1, max_nbases^2*femregion.nel]);
    jj_indexBd = reshape(vertcat(jj_indexBd{:}),[1, max_nbases^2*femregion.nel]);

    ii_indexVec = reshape(vertcat(ii_indexVec{:}),[1, max_nbases*femregion.nel]);

    Matrices = UpdateMatricesStructure(Matrices);

    % Extra-space in the matrices to be deleted
    select       = (jj_index > 0);
    selectBd     = (jj_indexBd > 0);
    selectVec    = (ii_indexVec > 0);
    select_neigh = (jj_index_neigh > 0);

    % Delete the extra-spaces in indices
    ii_index = ii_index(select);
    jj_index = jj_index(select);
    ii_indexBd = ii_indexBd(select);
    jj_indexBd = jj_indexBd(select);
    ii_index_neigh = ii_index_neigh(select_neigh);
    jj_index_neigh = jj_index_neigh(select_neigh);
    
    % Final reshape of the matrices
    if AssembInfo.assemblyvolume
        matFields = fieldnames(Matrices.Volume);
    % If we are assembling a multiphysics problem, apply the reshape to the substructs
 	if isstruct(Matrices.Volume.(matFields{1}))
            for kk = 1:length(matFields)
                Matrices.Volume.(matFields{kk}) = structfun(@(M) ReshapeVectorsMatrices(M, ii_index, jj_index, [], [], femregion, select, [], selectVec, AssembInfo.MatTag.(matFields{kk})), Matrices.Volume.(matFields{kk}), 'UniformOutput', false);
            end    
    % If we are assembling a single-physics problem, apply the reshape directly
        else
            Matrices.Volume = structfun(@(M) ReshapeVectorsMatrices(M, ii_index, jj_index, [], [], femregion, select, [], selectVec,[]), Matrices.Volume, 'UniformOutput', false);
        end
    end

    if AssembInfo.assemblyfaces
        matFields = fieldnames(Matrices.Faces);
    % If we are assembling a multiphysics problem, apply the reshape to the substructs
        if isstruct(Matrices.Faces.(matFields{1}))
            for kk = 1:length(matFields)
                Matrices.Faces.(matFields{kk}) = structfun(@(M) ReshapeVectorsMatrices(M, ii_indexBd, jj_indexBd, ii_index_neigh, jj_index_neigh, femregion, select, select_neigh, selectVec, AssembInfo.MatTag.(matFields{kk})), Matrices.Faces.(matFields{kk}), 'UniformOutput', false);
            end
    % If we are assembling a single-physics problem, apply the reshape directly
        else
            Matrices.Faces = structfun(@(M) ReshapeVectorsMatrices(M, ii_indexBd, jj_indexBd, ii_index_neigh, jj_index_neigh, femregion, select, select_neigh, selectVec, []), Matrices.Faces, 'UniformOutput', false);
        end
    end
    %% Final matrices for the solver
    [Matrices] = Funcs.FinalMatrices(Matrices);

end

