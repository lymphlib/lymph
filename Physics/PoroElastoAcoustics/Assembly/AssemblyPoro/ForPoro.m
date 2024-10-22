%> @file  ForPoro.m
%> @author Ilario Mazzieri
%> @date 25 June 2024
%> @brief Assembly of the rhs for the poroelastic problem \cite ABNM2021
%>
%==========================================================================
%> @section classForPoroPEA Class description
%==========================================================================
%> @brief Assembly of the rhs for the poroelastic problem \cite ABNM2021
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval F     Struct with the elsatic rhs terms
%>
%==========================================================================

function [F] = ForPoro(Data, neighbor, femregion)

%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the forcing term
F1 = zeros(femregion.ndof_p,1);
F2 = zeros(femregion.ndof_p,1);
J1 = zeros(femregion.ndof_p,1);
J2 = zeros(femregion.ndof_p,1);

G1 = zeros(femregion.ndof_p,1);
G2 = zeros(femregion.ndof_p,1);
H1 = zeros(femregion.ndof_p,1);
H2 = zeros(femregion.ndof_p,1);

% Dirichlet conditions
FD1 = zeros(femregion.ndof_p,1);
FD2 = zeros(femregion.ndof_p,1);
GD1 = zeros(femregion.ndof_p,1);
GD2 = zeros(femregion.ndof_p,1);



%% Loop over the elements

index_shift=0;

% Visualization of computational progress
prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);

for ie = 1:femregion.nel % loop over elements
    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);

    % Id and tag selection for elements
    tag_ie = femregion.tag(ie);

    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    index_element = index_shift + (1:1:femregion.nedges(ie))';
    index_shift   = index_element(end);

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

    % Check if the element elastic
    id_ie = femregion.id(ie);

    if tag_ie == 'P'

        % Local forcing vector definition
        F1_loc = zeros(femregion.nbases,1);
        F2_loc = zeros(femregion.nbases,1);
        J1_loc = zeros(femregion.nbases,1);
        J2_loc = zeros(femregion.nbases,1);
        G1_loc = zeros(femregion.nbases,1);
        G2_loc = zeros(femregion.nbases,1);
        H1_loc = zeros(femregion.nbases,1);
        H2_loc = zeros(femregion.nbases,1);

        % Local forcing term for Dirichelt conditions
        F1_diri_loc = zeros(femregion.nbases,1);
        F2_diri_loc = zeros(femregion.nbases,1);
        G1_diri_loc = zeros(femregion.nbases,1);
        G2_diri_loc = zeros(femregion.nbases,1);


        for iTria = 1:size(Tria,1)

            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Scaled weights
            dx = det(BJ) * w_2D;

            % Evaluation of forcing terms
            fSource1 = Data.source_up{1}(xq,yq);
            fSource2 = Data.source_up{2}(xq,yq);
            jSource1 = Data.source_upd{1}(xq,yq);
            jSource2 = Data.source_upd{2}(xq,yq);

            gSource1 = Data.source_wp{1}(xq,yq);
            gSource2 = Data.source_wp{2}(xq,yq);
            hSource1 = Data.source_wpd{1}(xq,yq);
            hSource2 = Data.source_wpd{2}(xq,yq);

            MxxSource = Data.sourceMxx_poro{1}(xq,yq);
            MxySource = Data.sourceMxy_poro{1}(xq,yq);
            MyxSource = Data.sourceMyx_poro{1}(xq,yq);
            MyySource = Data.sourceMyy_poro{1}(xq,yq);

            % Construction and evalutation on the quadrature points of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

            % Vector assembling
            F1_loc = F1_loc + (dx.*phiq)'*fSource1 + (dx.*gradqx)'*MxxSource + 0.5 * (dx.*gradqy)'*MxySource;
            F2_loc = F2_loc + (dx.*phiq)'*fSource2 + 0.5 * (dx.*gradqx)'*MyxSource + (dx.*gradqy)'*MyySource;
            J1_loc = J1_loc + (dx.*phiq)'*jSource1;
            J2_loc = J2_loc + (dx.*phiq)'*jSource2;

            G1_loc = G1_loc + (dx.*phiq)'*gSource1 + (dx.*gradqx)'*MxxSource + 0.5 * (dx.*gradqy)'*MxySource;
            G2_loc = G2_loc + (dx.*phiq)'*gSource2 + 0.5 * (dx.*gradqx)'*MyxSource + (dx.*gradqy)'*MyySource;
            H1_loc = H1_loc + (dx.*phiq)'*hSource1;
            H2_loc = H2_loc + (dx.*phiq)'*hSource2;

        end

        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

        % Loop over faces
        for iedg=1:neighbor.nedges(ie)

            % Dirichlet boundary faces
            if  neigh_ie(iedg) == -1

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
                mu  = Data.mu{id_ie}(xq,yq);
                lam = Data.lam{id_ie}(xq,yq);
                harm_ave = (lam+2*mu);
                beta = Data.beta{id_ie}(xq,yq);
                m = Data.m{id_ie}(xq,yq);

                % Auxiliary quantities (cf. physical parameters) for vector assembling
                aa = 0.5 * (lam+2*mu) * nx;
                ff = 0.5 * (lam+2*mu) * ny;
                bb = 0.5 * lam * nx;
                gg = 0.5 * lam * ny;
                ee = 0.5 * mu * nx;
                cc = 0.5 * mu * ny;
    
                % Construction and evalutation on the quadrature points of the basis functions
                [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

                % Evaluation of boundary functions
                gD1 = Data.DirBCPoro_up{1}(xq,yq);
                gD2 = Data.DirBCPoro_up{2}(xq,yq);
                hD1 = Data.DirBCPoro_wp{1}(xq,yq);
                hD2 = Data.DirBCPoro_wp{2}(xq,yq);

                % Vector assembling
                F1_diri_loc = F1_diri_loc - 2 * (ds .* (aa .* gradedgeqx + cc .* gradedgeqy))' *gD1 ...
                                - 2 * (ds .* (gg .* gradedgeqx + ee .* gradedgeqy))' *gD2 ;
                F1_diri_loc = F1_diri_loc + penalty_geom(iedg) * (ds .* harm_ave .* phiedgeq)' * gD1;

                F2_diri_loc = F2_diri_loc - 2 * (ds .* (cc .* gradedgeqx + bb .* gradedgeqy))' *gD1 ...
                                - 2 * (ds .* (ee .* gradedgeqx + ff .* gradedgeqy))' *gD2 ;
                F2_diri_loc = F2_diri_loc + penalty_geom(iedg) * (ds .* harm_ave .* phiedgeq)' * gD2;

                
                F1_diri_loc = F1_diri_loc - (ds .* (beta .* beta .*m .* gradedgeqx))' * (gD1*nx + gD2*ny);
                F1_diri_loc = F1_diri_loc + penalty_geom(iedg) * (ds .* (beta .* beta .*m .* nx .* nx .* phiedgeq ))' * gD1 ...
                                + penalty_geom(iedg) * (ds .* (beta .* beta .*m .* nx .* ny .* phiedgeq ))' * gD2;

                F1_diri_loc = F1_diri_loc - (ds .* (beta .*m .* gradedgeqx))' * (hD1*nx + hD2*ny);
                F1_diri_loc = F1_diri_loc + penalty_geom(iedg) * (ds .* (beta .*m .* nx .* nx .* phiedgeq ))' * hD1 ...
                                + penalty_geom(iedg) * (ds .* (beta .*m .* nx .* ny .* phiedgeq ))' * hD2;


                F2_diri_loc = F2_diri_loc - (ds .* (beta .* beta .*m .* gradedgeqy))' * (gD1*nx + gD2*ny);
                F2_diri_loc = F2_diri_loc + penalty_geom(iedg) * (ds .* (beta .* beta .*m .* nx .* ny .* phiedgeq ))' * gD1 ...
                                + penalty_geom(iedg) * (ds .* (beta .* beta .*m .* ny .* ny .* phiedgeq ))' * gD2;

                F2_diri_loc = F2_diri_loc - (ds .* (beta .*m .* gradedgeqy))' * (hD1*nx + hD2*ny);
                F2_diri_loc = F2_diri_loc + penalty_geom(iedg) * (ds .* (beta .*m .* nx .* ny .* phiedgeq ))' * hD1 ...
                                + penalty_geom(iedg) * (ds .* (beta .*m .* ny .* ny .* phiedgeq ))' * hD2;

 
                G1_diri_loc = G1_diri_loc - (ds .* (beta .*m .* gradedgeqx))' * (gD1*nx + gD2*ny);
                G1_diri_loc = G1_diri_loc + penalty_geom(iedg) * (ds .* (beta .*m .* nx .* nx .* phiedgeq ))' * gD1 ...
                                + penalty_geom(iedg) * (ds .* (beta .*m .* nx .* ny .* phiedgeq ))' * gD2;

                G1_diri_loc = G1_diri_loc - (ds .* (m .* gradedgeqx))' * (hD1*nx + hD2*ny);
                G1_diri_loc = G1_diri_loc + penalty_geom(iedg) * (ds .* (m .* nx .* nx .* phiedgeq ))' * hD1 ...
                                + penalty_geom(iedg) * (ds .* (m .* nx .* ny .* phiedgeq ))' * hD2;


                G2_diri_loc = G2_diri_loc - (ds .* (beta .*m .* gradedgeqy))' * (gD1*nx + gD2*ny);
                G2_diri_loc = G2_diri_loc + penalty_geom(iedg) * (ds .* (beta .*m .* nx .* ny .* phiedgeq ))' * gD1 ...
                                + penalty_geom(iedg) * (ds .* (beta .*m .* ny .* ny .* phiedgeq ))' * gD2;

                G2_diri_loc = G2_diri_loc - (ds .* (m .* gradedgeqy))' * (hD1*nx + hD2*ny);
                G2_diri_loc = G2_diri_loc + penalty_geom(iedg) * (ds .* (m .* nx .* ny .* phiedgeq ))' * hD1 ...
                                + penalty_geom(iedg) * (ds .* (m .* ny .* ny .* phiedgeq ))' * hD2;


            end
        end

        % Local vector to global vector
        F1(index) = F1_loc;
        F2(index) = F2_loc;
        J1(index) = J1_loc;
        J2(index) = J2_loc;
        G1(index) = G1_loc;
        G2(index) = G2_loc;
        H1(index) = H1_loc;
        H2(index) = H2_loc;

        FD1(index) = F1_diri_loc;
        FD2(index) = F2_diri_loc;
        GD1(index) = G1_diri_loc;
        GD2(index) = G2_diri_loc;

    end
end

F.f_p = [F1; F2];
F.g_p = [G1; G2];

F.j_p = [J1; J2];
F.h_p = [H1; H2];

F.f_p_diri = [FD1; FD2];
F.g_p_diri = [GD1; GD2];


fprintf('\n');

