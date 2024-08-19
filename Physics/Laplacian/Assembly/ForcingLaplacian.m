%> @file  ForcingLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Assembly of the RHS for the Poisson's problem
%>
%==========================================================================
%> @section classForcingLaplacian Class description
%==========================================================================
%> @brief            Assembly the RHS and boundary conditions for the
%>Poisson's problem
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval F         RHS vector
%>
%==========================================================================

function [F] = ForcingLaplacian(Data, neighbor, femregion)

    %% Quadrature values
    
    [ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);
    
    %% Initialization of the forcing term
    F = zeros(femregion.ndof,1);
    
    %% Loop over the elements
    
    % Visualization of computational progress
    prog = 0;
    fprintf(1,'Computation Progress: %3d%%\n',prog);
    
    for ie = 1:femregion.nel
        
        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);

        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        
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
    
        % Local forcing vector definition
        F_loc = zeros(femregion.nbases,1);

        for iTria = 1:size(Tria,1)
                 
             % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
        
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Scaled weights 
            dx = det(BJ) * w_2D;

            % Evaluation of physical parameters 
            fSource = Data.source{1}(xq,yq);
    
            % Construction of the basis functions
            phiq = Evalshape2D(femregion, ie, qNodes_2D);
            
            %% Vector assembling

            F_loc = F_loc + (dx.*phiq)'*fSource;
            
        end
            
        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);
        
        
        % Loop over faces
        for iedg = 1 : neighbor.nedges(ie)
               
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

            % Evaluation of physical parameters 
            mu = Data.mu{1}(xq,yq); 

            % Construction of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);
          
            % Extraction of normals to the face
            nx = normals(1,iedg);
            ny = normals(2,iedg);
              
            % Dirichlet boundary faces
            if neigh_ie(iedg) == -1
                            
               gD = Data.DirBC{1}(xq,yq);
     
               F_loc = F_loc - (ds .* mu .* ( nx * gradedgeqx + ny * gradedgeqy))' * gD;
               F_loc = F_loc + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * gD;
              
            end
        end

        % Local vector to global vector
        F(index) = F_loc;
    end
    
end
