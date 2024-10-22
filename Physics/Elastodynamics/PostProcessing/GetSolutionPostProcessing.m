%> @file  GetSolutionPostProcessing.m
%> @author Ilario Mazzieri, Stefano Bonetti
%> @date 9 October 2024
%> @brief  Compute solution for post-processing
%>
%==========================================================================
%> @section classElastodynamicsGetSolutionPostProcessing Class description
%==========================================================================
%> @brief  Compute solution for post-processing
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param Solution   Struct with solution vectors
%> @param time       time instant
%
%> @retval Xh        Cell array of length 1 containing a
%>                   struct with approximate solutions at quad nodes (order: uh, vh, ph)
%>                   and sequence of strings useful for post-processing
%> @retval Xexact    Cell array of length 1 containing a
%>                   struct with exact solutions at quad nodes (order: uex, vex, pex)
%>                   and sequence of strings useful for post-processing
%>
%==========================================================================

function [Xh, Xexact] = GetSolutionPostProcessing(Data, femregion, neighbor, Solution, time)

    %% Nodes values
    ref_qNodes_2D = linspace(0,1,Data.NPtsVisualization+2)';
    ref_qNodes_2D = ref_qNodes_2D(2:end-1);
    ref_qNodes_2D = [repmat(ref_qNodes_2D,Data.NPtsVisualization,1), reshape(repmat(ref_qNodes_2D,1,Data.NPtsVisualization)',Data.NPtsVisualization^2,1)];

%% Setup
    X  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Y  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Uh = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    Vh = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
    Ph = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Id = zeros(femregion.nel,2);
    Bd = cell(femregion.nel,1);

    if Data.PlotExact
        Uex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
        Vex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),2);
        Pex = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    end

    % shift index
    CountID = 0;
    id_shift_e = 0;

    %% Loop over the elements
    % Visualization of computational progress
    prog = 0;
    fprintf(1,'Computation Progress: %3d%%\n',prog);

    for ie = 1:femregion.nel % loop over elements

        % Visualization of computational progress
        prog = ( 100*(ie/femregion.nel) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);

        % Element id
        id_ie  = femregion.id(ie);

        % Selection of the matrix positions associated to element ie
        index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

        % Extraction of element geometrical information
        coords_ie = femregion.coords_element{ie};

        % Counter boundary nodes
        CountBD_min = 0;
        CountBD_max = 0;

        % First node of the current element
        Id(ie,1) = CountID + 1;

        % Elastic element
        if femregion.tag(ie) == 'E'

            index_e = index - femregion.ndof_p - femregion.ndof_a;

            % Local solution
            u1_loc = Solution(index_e);
            u2_loc = Solution(index_e + femregion.ndof_e);
            % Local velocity
            du1_loc = Solution(2*femregion.ndof_e + index_e);
            du2_loc = Solution(2*femregion.ndof_e + index_e + femregion.ndof_e);

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
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

            % Evaluation of physical parameters
            lambda = Data.lam_el{id_ie-id_shift_e}(xq,yq);

            % Approximated solutions at quadrature points
            ue1h_loc  = phiq*u1_loc;
            ue2h_loc  = phiq*u2_loc;
            peh_loc   = -lambda.*(gradqx*u1_loc + gradqy*u2_loc);

            due1h_loc  = phiq*du1_loc;
            due2h_loc  = phiq*du2_loc;

            % Exact and approximated solutions at quadrature points
            if Data.PlotExact
                ue1ex_loc  = Data.ue_ex{1}(xq,yq,time);
                ue2ex_loc  = Data.ue_ex{2}(xq,yq,time);
                due1ex_loc = Data.due_dt_ex{1}(xq,yq,time);
                due2ex_loc = Data.due_dt_ex{2}(xq,yq,time);
                pex_loc    = -lambda.*(Data.grad_ue_ex{1}(xq,yq,time) + Data.grad_ue_ex{4}(xq,yq,time));
            end

            % Fill the output (Displacement, Pressure, and Velocity)
            X(CountID + 1 : CountID + lqn,:)   = xq;
            Y(CountID + 1 : CountID + lqn,:)   = yq;

            Uh(CountID + 1 : CountID + lqn,:)  = [ue1h_loc, ue2h_loc];
            Ph(CountID + 1 : CountID + lqn,:)  = peh_loc;
            Vh(CountID + 1 : CountID + lqn,:)  = [due1h_loc, due2h_loc];

            if Data.PlotExact
                Uex(CountID + 1 : CountID + lqn,:) = [ue1ex_loc, ue2ex_loc];
                Pex(CountID + 1 : CountID + lqn,:) = pex_loc;
                Vex(CountID + 1 : CountID + lqn,:) = [due1ex_loc, due2ex_loc];
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

                % Evaluation of physical parameters
                lambda = Data.lam_el{id_ie-id_shift_e}(xq,yq);

                % Construction and evaluation of the basis functions in the visualization points
                [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

                % Approximated solutions at quadrature points
                ue1h_loc  = phiedgeq*u1_loc;
                ue2h_loc  = phiedgeq*u2_loc;
                peh_loc   = -lambda.*(gradedgeqx*u1_loc + gradedgeqy*u2_loc);

                due1h_loc  = phiedgeq*du1_loc;
                due2h_loc  = phiedgeq*du2_loc;

                % Fill the output
                X(CountID + 1 : CountID + lqn,:)   = xq;
                Y(CountID + 1 : CountID + lqn,:)   = yq;
                Uh(CountID + 1 : CountID + lqn,:)  = [ue1h_loc, ue2h_loc];
                Ph(CountID + 1 : CountID + lqn,:)  = peh_loc;
                Vh(CountID + 1 : CountID + lqn,:)  = [due1h_loc, due2h_loc];

                if Data.PlotExact
                    ue1ex_loc  = Data.ue_ex{1}(xq,yq,time);
                    ue2ex_loc  = Data.ue_ex{2}(xq,yq,time);
                    due1ex_loc = Data.due_dt_ex{1}(xq,yq,time);
                    due2ex_loc = Data.due_dt_ex{2}(xq,yq,time);
                    pex_loc    = -lambda.*(Data.grad_ue_ex{1}(xq,yq,time) + Data.grad_ue_ex{4}(xq,yq,time));

                    % Fill the output
                    Uex(CountID + 1 : CountID + lqn,:) = [ue1ex_loc, ue2ex_loc];
                    Pex(CountID + 1 : CountID + lqn,:) = pex_loc;
                    Vex(CountID + 1 : CountID + lqn,:) = [due1ex_loc, due2ex_loc];
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

    XPlot{1} = X(1:CountID);
    XPlot{2} = Y(1:CountID);
    XPlot{3} = Uh(1:CountID,:);
    XPlot{4} = Vh(1:CountID,:);
    XPlot{5} = Ph(1:CountID,:);
    StrVTK   = {'x','y','uh','vh','ph'};
    StrPlot  = {'x','y','u','v','p'};
    StrCSV   = {'x','y','uh1','uh2','vh1','vh2','ph'};

    Xh = {};

    Xh{1}.Solution = XPlot;
    Xh{1}.StrVTK   = StrVTK;
    Xh{1}.StrPlot  = StrPlot;
    Xh{1}.StrCSV   = StrCSV;
    Xh{1}.Id       = Id;
    Xh{1}.Bd       = Bd;

    if Data.PlotExact
        XPlotExact{1} = X(1:CountID);
        XPlotExact{2} = Y(1:CountID);
        XPlotExact{3} = Uex(1:CountID,:);
        XPlotExact{4} = Vex(1:CountID,:);
        XPlotExact{5} = Pex(1:CountID,:);

        Xexact = {};
        Xexact{1}.Solution = XPlotExact;
        Xexact{1}.StrPlot  = StrPlot;
    else
        Xexact = {[]};
    end

end
