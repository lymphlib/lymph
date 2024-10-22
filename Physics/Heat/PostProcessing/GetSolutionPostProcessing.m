%> @file  GetSolutionPostProcessing.m
%> @author Mattia Corti, Stefano Bonetti
%> @date 9 October 2024
%> @brief Compute solution for post-processing
%>
%==========================================================================
%> @section classHeatGetSolutionPostProcessing Class description
%==========================================================================
%> @brief Compute solution for post-processing
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param u_h        Numerical solution to evaluate in quadrature nodes
%> @param t          Time associated to the forcing term
%>
%> @retval Xh        Cell array of length 1 containing a
%>                   struct with approximate solutions at quad nodes and
%>                   sequence of strings useful for post-processing
%> @retval Xexact    Cell array of length 1 containing a
%>                   struct with exact solutions at quad nodes and 
%>                   sequence of strings useful for post-processing
%>
%==========================================================================

function [Xh, Xexact] = GetSolutionPostProcessing(Data, femregion, neighbor, u_h, t)
    
    %% Nodes values
    ref_qNodes_2D = linspace(0,1,Data.NPtsVisualization+2)';
    ref_qNodes_2D = ref_qNodes_2D(2:end-1);
    ref_qNodes_2D = [repmat(ref_qNodes_2D,Data.NPtsVisualization,1), reshape(repmat(ref_qNodes_2D,1,Data.NPtsVisualization)',Data.NPtsVisualization^2,1)];

    %% Setup
    X  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Y  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Uh = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Id = zeros(femregion.nel,2);
    Bd = cell(femregion.nel,1);
    
    if Data.PlotExact
        Uex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    end
     
    % shift index
    CountID = 0;

    %% Loop over the elements
    % Visualization of computational progress
    prog = 0;
    fprintf(1,'Computation Progress: %3d%%\n',prog);
    
    for ie = 1:femregion.nel % loop over elements
        
        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        
        % Extraction of element geometrical information
        coords_ie = femregion.coords_element{ie};

        % Counter boundary nodes
        CountBD_min = 0;
        CountBD_max = 0;
        
        % First node of the current element
        Id(ie,1) = CountID + 1;
        
        % Computation of nodes on the physical element
        qNodes_2D = ref_qNodes_2D;

        qNodes_2D(:,1) = ref_qNodes_2D(:,1)*(femregion.bbox(ie,2)-femregion.bbox(ie,1))+femregion.bbox(ie,1);
        qNodes_2D(:,2) = ref_qNodes_2D(:,2)*(femregion.bbox(ie,4)-femregion.bbox(ie,3))+femregion.bbox(ie,3);
        
        if Data.NPtsVisualization > 1
            dist_2D = norm(qNodes_2D(end,:)-qNodes_2D(end-1,:));
        end

        qNodes_2D = qNodes_2D(inpolygon(qNodes_2D(:,1),qNodes_2D(:,2),coords_ie(:,1),coords_ie(:,2)),:);

        xq  = qNodes_2D(:,1);
        yq  = qNodes_2D(:,2);
        lqn = length(xq);

        % Construction and evaluation of the basis functions in the visualization points
        phiq = Evalshape2D(femregion, ie, qNodes_2D);
            
        % Approximated solutions at quadrature points
        uh_loc = phiq*u_h(index);

        % Fill the output
        X(CountID+1:CountID+lqn, :) = xq;
        Y(CountID+1:CountID+lqn, :) = yq;
        Uh(CountID+1:CountID+lqn,:) = uh_loc;

        if Data.PlotExact
            % Evaluation of exact solution
            local_exact_u = Data.u_ex(xq,yq,t);

            % Fill the output
            Uex(CountID+1:CountID+lqn,:) = local_exact_u;
        end

        % Update auxiliary indexes
        CountID = CountID + lqn;
        CountBD_min = CountBD_min + lqn + 1;
        CountBD_max = CountBD_max + lqn;

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

            if Data.NPtsVisualization > 1
                ref_qNodes_1D = linspace(0,1,ceil(norm(p1-p2)/dist_2D)+1)';
                ref_qNodes_1D = ref_qNodes_1D(1:end-1);

                [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);
            else
                qNodes_1D = p1;
            end
 
            % Construction of quadrature nodes on the face
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            lqn = length(xq);

            % Construction and evaluation of the basis functions in the visualization points
            phiedgeq = Evalshape2D(femregion, ie, qNodes_1D);

            % Approximated solutions at quadrature points
            uh_loc = phiedgeq*u_h(index);
            
            % Fill the output
            X(CountID+1:CountID+lqn, :) = xq;
            Y(CountID+1:CountID+lqn, :) = yq;
            Uh(CountID+1:CountID+lqn,:) = uh_loc;
            
            if Data.PlotExact
                % Evaluation of exact solution
                local_exact_u = Data.u_ex(xq,yq,t);
                
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

    Xh = {};

    Xh{1}.StrVTK   = {'x', 'y', 'uh'};
    Xh{1}.StrCSV   = {'x', 'y', 'uh'};
    Xh{1}.StrPlot  = {'x', 'y', 'u'};
    Xh{1}.Solution = { X(1:CountID),   Y(1:CountID),   Uh(1:CountID)};
    Xh{1}.Id       = Id;
    Xh{1}.Bd       = Bd;

    if Data.PlotExact
        Xexact = {};
        Xexact{1}.StrPlot  = {'x', 'y', 'u'};
        Xexact{1}.Solution = { X(1:CountID),   Y(1:CountID),   Uex(1:CountID)};
    else
        Xexact = {[]};
    end

end
