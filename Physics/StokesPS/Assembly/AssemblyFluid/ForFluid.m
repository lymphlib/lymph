%> @file  ForFluid.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Assembly of the rhs for the Stokes problem
%>
%==========================================================================
%> @section classForFluid Class description
%==========================================================================
%> @brief Assembly of the rhs for the Stokes problem
%>
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%>
%> @retval F     Struct with the rhs terms
%>
%==========================================================================

function [F] = ForFluid(Data, neighbor, femregion)

%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the forcing term

f11 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
f12 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
f21 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
f22 = spalloc(femregion.ndof_f,1,femregion.ndof_f);

g11 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
g12 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
g21 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
g22 = spalloc(femregion.ndof_f,1,femregion.ndof_f);

%% Loop over the elements

index_shift = 0;

% Visualization of computational progress
prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);


for ie = 1:femregion.nel % loop over elements
    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    % Id and tag selection for elements
    id_ie  = femregion.id(ie);
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
    
    % Check if the element is fluid
    if tag_ie == 'F'
        
        index_f = index - femregion.ndof_p;
        
        % Local forcing vector definition
        f11_loc = zeros(femregion.nbases,1);
        f12_loc = zeros(femregion.nbases,1);
        f21_loc = zeros(femregion.nbases,1);
        f22_loc = zeros(femregion.nbases,1);
        g11_loc = zeros(femregion.nbases,1);
        g12_loc = zeros(femregion.nbases,1);
        g21_loc = zeros(femregion.nbases,1);
        g22_loc = zeros(femregion.nbases,1);
        
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Evaluation of forcing terms
            f11Source = Data.source_sigma{1}(xq,yq);
            f12Source = Data.source_sigma{2}(xq,yq);
            f21Source = Data.source_sigma{3}(xq,yq);
            f22Source = Data.source_sigma{4}(xq,yq);
            
            g11Source = Data.source_sigma_d{1}(xq,yq);
            g12Source = Data.source_sigma_d{2}(xq,yq);
            g21Source = Data.source_sigma_d{3}(xq,yq);
            g22Source = Data.source_sigma_d{4}(xq,yq);
            
            % Construction of the basis functions
            phiq = Evalshape2D(femregion, ie, qNodes_2D);
            
            %% Vector assembling
            
            f11_loc = f11_loc + (dx.*phiq)'*f11Source;
            f12_loc = f12_loc + (dx.*phiq)'*f12Source;
            f21_loc = f21_loc + (dx.*phiq)'*f21Source;
            f22_loc = f22_loc + (dx.*phiq)'*f22Source;
            
            g11_loc = g11_loc + (dx.*phiq)'*g11Source;
            g12_loc = g12_loc + (dx.*phiq)'*g12Source;
            g21_loc = g21_loc + (dx.*phiq)'*g21Source;
            g22_loc = g22_loc + (dx.*phiq)'*g22Source;
            
        end
        
        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);
        
        
        % Loop over faces
        for iedg=1:neighbor.nedges(ie) % loop over faces
            
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
            
            % Construction of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);
            
            % Neumann boundary condition for sigma -> div(sigma) = gd
            if  neigh_ie(iedg) == -1
                
                
                gD1 = Data.NeuBCsigma_x{1}(xq,yq) + Data.NeuBCsigma_y{2}(xq,yq);
                gD2 = Data.NeuBCsigma_x{3}(xq,yq) + Data.NeuBCsigma_y{4}(xq,yq);
                
                f11_loc = f11_loc + (ds.* phiedgeq)' * gD1 * nx;
                f12_loc = f12_loc + (ds.* phiedgeq)' * gD1 * ny;
                f21_loc = f21_loc + (ds.* phiedgeq)' * gD2 * nx;
                f22_loc = f22_loc + (ds.* phiedgeq)' * gD2 * ny;
                
            % Dirichlet boundary condition for sigma -> sigma.n = gn
            elseif  neigh_ie(iedg) == -2
                
                gN11 = Data.DirBCsigma{1}(xq,yq) * nx + Data.DirBCsigma{2}(xq,yq) * ny;
                gN12 = Data.DirBCsigma{1}(xq,yq) * nx + Data.DirBCsigma{2}(xq,yq) * ny;
                gN21 = Data.DirBCsigma{3}(xq,yq) * nx + Data.DirBCsigma{4}(xq,yq) * ny;
                gN22 = Data.DirBCsigma{3}(xq,yq) * nx + Data.DirBCsigma{4}(xq,yq) * ny;
                
                
                f11_loc = f11_loc + penalty_geom(iedg) * (ds.* phiedgeq)' * gN11 * nx - (ds.* gradedgeqx)' * gN11;
                f12_loc = f12_loc + penalty_geom(iedg) * (ds.* phiedgeq)' * gN12 * ny - (ds.* gradedgeqy)' * gN12;
                f21_loc = f21_loc + penalty_geom(iedg) * (ds.* phiedgeq)' * gN21 * nx - (ds.* gradedgeqx)' * gN21;
                f22_loc = f22_loc + penalty_geom(iedg) * (ds.* phiedgeq)' * gN22 * ny - (ds.* gradedgeqy)' * gN22;
                
                
            end
            
        end
        
        % Local vector to global vector
        f11(index_f) = f11_loc;
        f12(index_f) = f12_loc;
        f21(index_f) = f21_loc;
        f22(index_f) = f22_loc;
        g11(index_f) = g11_loc;
        g12(index_f) = g12_loc;
        g21(index_f) = g21_loc;
        g22(index_f) = g22_loc;

    end
end

F.f_f = [f11; f12; f21; f22];
F.g_f = [g11; g12; g21; g22];

fprintf('\n');


