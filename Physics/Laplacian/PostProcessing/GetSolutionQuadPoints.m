%> @file  GetSolutionQuadPoints.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti
%> @date 1 August 2024
%> @brief Compute the solution at quadrature nodes.
%>
%==========================================================================
%> @section classGetSolutionQuadPoints Class description
%==========================================================================
%> @brief  Compute the solution at quadrature nodes.
%>
%> @param Data        Struct with problem's data
%> @param femregion   Struct for finite elements data
%> @param neighbor    Struct for neighboring elements data
%> @param U           Modal solution
%>
%> @retval Xh         Cell array of length 1 containing a
%>                     struct with approximate solutions at quad nodes and
%>                     sequence of strings useful for post-processing
%> @retval Xexact     Cell array of length 1 containing a
%>                     struct with exact solutions at quad nodes and 
%>                     sequence of strings useful for post-processing
%>
%==========================================================================

function [Xh, Xexact] = GetSolutionQuadPoints(Data, femregion, neighbor, U)

%% Quadrature values
[ref_qNodes_1D, ~, ref_qNodes_2D, ~] = Quadrature(Data.NqnVisualization);

%% Setup
X  = zeros(Data.NqnVisualization.^2*sum(neighbor.nedges-2)+sum(neighbor.nedges),1);
Y  = zeros(Data.NqnVisualization.^2*sum(neighbor.nedges-2)+sum(neighbor.nedges),1);
Uh = zeros(Data.NqnVisualization.^2*sum(neighbor.nedges-2)+sum(neighbor.nedges),1);
Id = zeros(femregion.nel,2);
Bd = cell(femregion.nel,1);
    
if Data.PlotExact
    Uex = zeros(Data.NqnVisualization.^2*sum(neighbor.nedges-2)+sum(neighbor.nedges),1);
end
   
% shift index
CountID = 0;

%% Loop over the elements
% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie=1:femregion.nel
        
        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        
        % Extraction of element geometrical information
        coords_ie = femregion.coords_element{ie};

        % Creation of the subtriangulation of the element
        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

        % Counter boundary nodes
        CountBD_min = 0;
        CountBD_max = 0;
        
        % First node of the current element
        Id(ie,1) = CountID + 1;
    
    if femregion.tag(ie) == 'L'
        % Local solution
        u_loc = U(index);
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [~, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            lqn = length(xq);
         
            % Construction of the basis functions
            [phiq, ~, ~] = Evalshape2D(femregion, ie, qNodes_2D);
            
            % Approximated solutions at quadrature points
            uh_loc  = phiq*u_loc;
            
            % Fill the output
            X(CountID + 1 : CountID + lqn, :) = xq;
            Y(CountID + 1 : CountID + lqn, :) = yq;
            Uh(CountID + 1 : CountID + lqn,:) = uh_loc;

            if Data.PlotExact
                % Evaluation of exact solution
                uex_loc = Data.u_ex{1}(xq,yq);
      
                % Fill the output
                Uex(CountID + 1 : CountID + lqn, :) = uex_loc;
            end

            % Update auxiliary indexes
            CountID = CountID + lqn;
            CountBD_min = CountBD_min + lqn;
            CountBD_max = CountBD_max + lqn; 
            
        end

        % First boundary node of the current element
        CountBD_min = CountBD_min+1;
                
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

            [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], [0; ref_qNodes_1D]);
            
            % Construction of quadrature nodes on the face
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            lqn = length(xq);

            % Construction of the basis functions
            phiedgeq = Evalshape2D(femregion, ie, qNodes_1D);

            % Approximated solutions at quadrature points
            uh_loc = phiedgeq*u_loc;
          
            % Fill the output
            X(CountID+1:CountID+lqn, :) = xq;
            Y(CountID+1:CountID+lqn, :) = yq;
            Uh(CountID+1:CountID+lqn,:) = uh_loc;
            
            if Data.PlotExact
                % Evaluation of exact solution
                local_exact_u = Data.u_ex{1}(xq,yq);
                
                % Fill the output
                Uex(CountID+1:CountID+lqn,:) = local_exact_u;
            end

            % Update auxiliary indexes
            CountID  = CountID  + lqn;
            CountBD_max = CountBD_max + lqn; 
        end

        % save a vector from first to last (index of) boundary node of the current element
        Bd{ie} = (CountBD_min:CountBD_max)';

        % Last node of the current element
        Id(ie,2) = CountID;

    end
end

Xh = {};

Xh{1}.StrVTK   = {'x', 'y', 'uh'};
Xh{1}.StrCSV   = {'x', 'y', 'uh'};
Xh{1}.StrPlot  = {'x', 'y', 'u'};
Xh{1}.Solution = { X,   Y,   Uh};
Xh{1}.Id       = Id;
Xh{1}.Bd       = Bd;

if Data.PlotExact
    Xexact = {};

    Xexact{1}.StrPlot  = {'x', 'y', 'u'};
    Xexact{1}.Solution = { X,   Y,   Uex};
else
    Xexact = {[]};
end
