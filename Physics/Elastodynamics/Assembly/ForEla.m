%> @file  ForEla.m
%> @author Ilario Mazzieri, Stefano Bonetti
%> @date 24 July 2024
%> @brief Assembly of the rhs for the elastic problem \cite AM2017
%>
%==========================================================================
%> @section classForEla Class description
%==========================================================================
%> @brief Assembly of the rhs for the elastic problem \cite AM2017
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param time       Time instant
%
%> @retval F     Struct with the elsatic rhs terms
%>
%==========================================================================

function [F, G] = ForEla(Data, neighbor, femregion, time)

%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the forcing term
F1 = zeros(femregion.ndof_e,1);
F2 = zeros(femregion.ndof_e,1);
G1 = zeros(femregion.ndof_e,1);
G2 = zeros(femregion.ndof_e,1);

%% Loop over the elements

index_shift=0;
id_shift = 0;

% Visualization of computational progress
prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);

for ie=1:femregion.nel % loop over elements
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
    
    if tag_ie == 'E'
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        
        % Local forcing vector definition
        F1_loc = zeros(femregion.nbases,1);
        F2_loc = zeros(femregion.nbases,1);
        G1_loc = zeros(femregion.nbases,1);
        G2_loc = zeros(femregion.nbases,1);
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Evaluation of forcing terms
            fSource1 = Data.source_f{1}(xq,yq,time);
            fSource2 = Data.source_f{2}(xq,yq,time);
            gSource1 = Data.source_g{1}(xq,yq,time);
            gSource2 = Data.source_g{2}(xq,yq,time);
            
            MxxSource = Data.sourceMxx_el{1}(xq,yq,time);
            MxySource = Data.sourceMxy_el{1}(xq,yq,time);
            MyxSource = Data.sourceMyx_el{1}(xq,yq,time);
            MyySource = Data.sourceMyy_el{1}(xq,yq,time);
            
            % Construction and evalutation on the quadrature points of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
                        
            % Vector assembling
            F1_loc = F1_loc + (dx.*phiq)'*fSource1 + (dx.*gradqx)'*MxxSource + 0.5 * (dx.*gradqy)'*MxySource;
            F2_loc = F2_loc + (dx.*phiq)'*fSource2 + 0.5 * (dx.*gradqx)'*MyxSource + (dx.*gradqy)'*MyySource;
            G1_loc = G1_loc + (dx.*phiq)'*gSource1;
            G2_loc = G2_loc + (dx.*phiq)'*gSource2;
            
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
                mu  = Data.mu_el{id_ie-id_shift}(xq,yq);
                lam = Data.lam_el{id_ie-id_shift}(xq,yq);
                harm_ave = (lam+2*mu);
                
                % Auxiliary quantities (cf. physical parameters) for vector assembling
                aa = (lam+2*mu) * nx;
                ff = (lam+2*mu) * ny;
                bb = lam * nx;
                gg = lam * ny;
                ee = mu * nx;
                cc = mu * ny;
                
                % Construction and evalutation on the quadrature points of the basis functions
                [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);
                
                % Evaluation of boundary functions
                gD1 = Data.DirBCEla{1}(xq,yq,time);
                gD2 = Data.DirBCEla{2}(xq,yq,time);
                
                % Vector assembling
                F1_loc = F1_loc - (ds.* ( aa .* gradedgeqx + cc .* gradedgeqy))' * gD1 - (ds.* ( gg .* gradedgeqx + ee .* gradedgeqy))' * gD2;
                F1_loc = F1_loc + penalty_geom(iedg) * (ds .* harm_ave .* phiedgeq)' * gD1;
                
                F2_loc = F2_loc - (ds.* ( cc .* gradedgeqx + bb .* gradedgeqy))' * gD1 - (ds.* ( ee .* gradedgeqx + ff .* gradedgeqy))' * gD2;
                F2_loc = F2_loc + penalty_geom(iedg) * (ds .* harm_ave .* phiedgeq)' * gD2;
                
            end
        end
        
        % Local vector to global vector
        F1(index_e) = F1_loc;
        F2(index_e) = F2_loc;
        G1(index_e) = G1_loc;
        G2(index_e) = G2_loc;
               
    end
end

F = [F1; F2];
G = [G1; G2];

fprintf('\n');

