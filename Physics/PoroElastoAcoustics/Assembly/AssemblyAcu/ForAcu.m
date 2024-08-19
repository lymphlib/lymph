%> @file  ForAcu.m
%> @author Ilario Mazzieri
%> @date 25 June 2024
%> @brief Assembly of the rhs for the acoustic problem \cite ABNM2021
%>
%==========================================================================
%> @section classForAcu Class description
%==========================================================================
%> @brief Assembly of the rhs for the acoustic problem \cite ABNM2021
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param F          Struct with rhs terms
%>
%> @retval F     Struct with rhs terms
%>
%==========================================================================

function [F] = ForAcu(Data, neighbor, femregion, F)

%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the forcing term
F1 = zeros(femregion.ndof_a,1);
FD1 = zeros(femregion.ndof_a,1);

%% Loop over the elements

index_shift=0;
id_shift = max(Data.TagElPoro);
if(isempty(id_shift)); id_shift = 0; end

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
    
    if tag_ie == 'A'
        
        index_a = index - femregion.ndof_p;
        
        % Local forcing vector definition
        F1_loc = zeros(femregion.nbases,1);
        F1_diri_loc = zeros(femregion.nbases,1);
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Evaluation of forcing terms
            rho_a = Data.rho_a{id_ie-id_shift}(xq,yq);
            fSource1 = Data.source_phi{1}(xq,yq);
            
            % Construction of the basis functions
            phiq = Evalshape2D(femregion, ie, qNodes_2D);
            
            
            %% Vector assembling
            F1_loc = F1_loc + (dx .* rho_a .* phiq)' * fSource1;
            
        end
        
        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);
        
        % Loop over faces
        for iedg = 1:neighbor.nedges(ie)
            
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
            rho_a = Data.rho_a{id_ie-id_shift}(xq,yq);
            
            aa = rho_a * nx;
            bb = rho_a * ny;
            
            % Construction of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);
            
            
            % Dirichlet boundary faces
            if  neigh_ie(iedg) == -1
                
                gD1 = Data.DirBCAcu{1}(xq,yq);
                
                F1_diri_loc = F1_diri_loc - ( ds.* ( aa .* gradedgeqx + bb .* gradedgeqy))' * gD1;
                F1_diri_loc = F1_diri_loc + penalty_geom(iedg) * (ds .* rho_a .* phiedgeq)' * gD1;
                
                
            end
        end
        
        % Local vector to global vector
        F1(index_a) = F1_loc;
        FD1(index_a) = F1_diri_loc;
        
        
    end
end

F.f_a = F1;
F.f_a_diri = FD1;

fprintf('\n');

