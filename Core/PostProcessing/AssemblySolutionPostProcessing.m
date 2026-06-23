%> @file  AssemblySolutionPostProcessing.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti, Caterina Leimer Saglio
%> @date 5 June 2026
%> @brief Compute the solution for post-processing.
%>
%==========================================================================
%> @section AssemblySolutionPostProcessing Class description
%==========================================================================
%> @brief  Compute the solution for post-processing of single- and multi-physics problems.
%>
%> @param Data          Struct with problem's data
%> @param femregion     Struct for finite elements data
%> @param neighbor      Struct for neighboring elements data
%> @param SolutionInfo  Cell array containing quantity to postprocess and
%> the necessary information. For each cell element, it must contain:
%> - outputnames: names of the solution/properties (i.e. degree) to export
%> - outputsizes: 1 if it is a scalar or 2 if it is a vector
%> - outputsol: vector containing cell values or modal coefficients
%> - outputexpr: expression to recover the desired output from the
%> outputsol and the basis functions and gradients (phiq, gradphiqx,
%> gradphiqy)
%> - outputexact: expression of the exact value associated to the quantity 
%> in terms of coordinates (xq,yq) and time (t)
%> - t: time variable (= 0 if it is a steady problem)
%>
%> @retval Xh         Cell array containing a struct with approximate
%> solutions and strings useful for post processing. Each cell of the array 
%> contains a different physics for the multiphysics implementation 
%> compatibility. The struct contains:
%> - StrVTK: name of the variables in the VTK output
%> - StrPlot: name of the variables in the scatter plot output   
%> - StrCSV: name of the variables in the CSV output
%> - Solution: Vectors containing the coordinates (x,y) and the solution values
%> - Id: Matrix containing for the j-th mesh element the first and last 
%>   indexes of the coordinates and solution vectors
%> - Bd: Cell array containing the id of the coordinates on the cell boundary
%> @retval Xexact     Cell array containing a struct with exact
%> solutions and strings useful for post processing. Each cell of the array 
%> contains a different physics for the multiphysics implementation (see Xh
%> for detailed structure)
%>
%==========================================================================

function [Xh, Xexact] = AssemblySolutionPostProcessing(Data, femregion, neighbor, SolutionInfo)

    %% Nodes values
    ref_qNodes_2D = linspace(0,1,Data.NPtsVisualization+2)';
    ref_qNodes_2D = ref_qNodes_2D(2:end-1);
    ref_qNodes_2D = [repmat(ref_qNodes_2D,Data.NPtsVisualization,1), reshape(repmat(ref_qNodes_2D,1,Data.NPtsVisualization)',Data.NPtsVisualization^2,1)];

    %% Setup
    X  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    Y  = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),1);
    
    Solh = cell(size(SolutionInfo.outputnames));
    
    for ii = 1:length(SolutionInfo.outputnames)
        Solh{ii} = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),SolutionInfo.outputsizes{ii});
    end

    Id = zeros(femregion.nel,2);
    Bd = cell(femregion.nel,1);

    if Data.PlotExact
        Solex = cell(length(SolutionInfo.outputnames));

        for ii = 1:length(SolutionInfo.outputnames)
            Solex{ii} = zeros(Data.NPtsVisualization.^2*femregion.nel+(Data.NPtsVisualization+1)*sum(neighbor.nedges),SolutionInfo.outputsizes{ii});
        end
    end

    % shift index
    CountID = 0;

    %% Loop over the elements
    for ie=1:femregion.nel

        if femregion.label(ie) == SolutionInfo.label

            % Selection of the matrix positions associated to element ie
            index = (sum(femregion.nbases(1:ie-1)))* ones(femregion.nbases(ie),1) + (1:femregion.nbases(ie))';

            % Extraction of element geometrical information
            coords_ie = femregion.coords_element{ie};

            % Counter boundary nodes
            CountBD_min = 0;
            CountBD_max = 0;

            % First node of the current element
            Id(ie,1) = CountID + 1;

            % Computation of 2D-nodes on the physical element
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

            % Fill the output
            X(CountID + 1 : CountID + lqn, :) = xq;
            Y(CountID + 1 : CountID + lqn, :) = yq;

            % Evaluate eventual physical parameters
            if isfield(SolutionInfo,'phys_param')
                phys_param = zeros(lqn,length(SolutionInfo.phys_param));
                for ii=1:length(SolutionInfo.phys_param)
                    phys_param(1:lqn,ii) = SolutionInfo.phys_param{ii}(xq,yq,SolutionInfo.t,femregion.id_phys(ie));
                end
            end

            for ii = 1:length(SolutionInfo.outputnames)
                if length(SolutionInfo.outputsol{ii}) == femregion.nel
                    eval_idx = ie;
                    if isfield(SolutionInfo,'tag')
                        eval_idx = eval_idx - sum(femregion.nel_phys(1:find(femregion.phys == SolutionInfo.label)-1));
                    end
                else
                    eval_idx = index;
                    if isfield(SolutionInfo,'tag')
                        eval_idx = eval_idx - sum(femregion.ndof_phys(1:find(femregion.phys == SolutionInfo.label)-1));
                    end
                end

                SolOut = SolutionInfo.outputsol{ii}(eval_idx,:);

                if nargin(SolutionInfo.outputexpr{ii}) == 4
                    Solh{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexpr{ii}(phiq, gradqx, gradqy, SolOut);
                else
                    Solh{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexpr{ii}(phiq, gradqx, gradqy, SolOut, phys_param);
                end

                if Data.PlotExact
                    if nargin(SolutionInfo.outputexact{ii}) == 3
                        Solex{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexact{ii}(xq,yq,SolutionInfo.t);
                    else
                        Solex{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexact{ii}(xq,yq,SolutionInfo.t, phys_param);
                    end
                end
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

                % Computation of 1D-nodes on the physical element
                if Data.NPtsVisualization > 1
                    ref_qNodes_1D = linspace(0,1,ceil(norm(p1-p2)/dist_2D)+1)';
                    ref_qNodes_1D = ref_qNodes_1D(1:end-1);

                    [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);
                else
                    qNodes_1D = p1;
                end

                % Construction of evaluation nodes on the face
                xq  = qNodes_1D(:,1);
                yq  = qNodes_1D(:,2);
                lqn = length(xq);

                % Construction and evaluation of the basis functions in the visualization points
                [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

                % Fill the output
                X(CountID+1:CountID+lqn, :) = xq;
                Y(CountID+1:CountID+lqn, :) = yq;

                % Evaluate eventual physical parameters
                if isfield(SolutionInfo,'phys_param')
                    phys_param = zeros(lqn,length(SolutionInfo.phys_param));
                    for ii=1:length(SolutionInfo.phys_param)
                        phys_param(1:lqn,ii) = SolutionInfo.phys_param{ii}(xq,yq,SolutionInfo.t,femregion.id_phys(ie));
                    end
                end

                for ii = 1:length(SolutionInfo.outputnames)
                    if length(SolutionInfo.outputsol{ii}) == femregion.nel
                        eval_idx = ie;
                        if isfield(SolutionInfo,'tag')
                            eval_idx = eval_idx - sum(femregion.nel_phys(1:find(femregion.phys == SolutionInfo.label)-1));
                        end
                    else
                        eval_idx = index;
                        if isfield(SolutionInfo,'tag')
                            eval_idx = eval_idx - sum(femregion.ndof_phys(1:find(femregion.phys == SolutionInfo.label)-1));
                        end                   
                    end

                    SolOut = SolutionInfo.outputsol{ii}(eval_idx,:);

                    if nargin(SolutionInfo.outputexpr{ii}) == 4
                        Solh{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexpr{ii}(phiedgeq, gradedgeqx, gradedgeqy, SolOut);
                    else
                        Solh{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexpr{ii}(phiedgeq, gradedgeqx, gradedgeqy, SolOut, phys_param);
                    end

                    if Data.PlotExact
                        if nargin(SolutionInfo.outputexact{ii}) == 3
                            Solex{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexact{ii}(xq,yq,SolutionInfo.t);
                        else
                            Solex{ii}(CountID + 1 : CountID + lqn,:) = SolutionInfo.outputexact{ii}(xq,yq,SolutionInfo.t, phys_param);
                        end
                    end
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

    %% Construction of the output solution structure
    Xh = {};

    Xh{1}.StrVTK   = {'x', 'y'};
    Xh{1}.StrPlot  = {'x', 'y'};
    Xh{1}.StrCSV   = {'x', 'y'};

    Xh{1}.Solution = { X(1:CountID),   Y(1:CountID)};
    Xh{1}.Id       = Id(all(Id>0,2),:);
    Xh{1}.Bd       = Bd(~cellfun('isempty', Bd));

    for ii = 1:length(SolutionInfo.outputnames)

        % Construction of names to save the VTK file
        Xh{1}.StrVTK{end+1}  = [SolutionInfo.outputnames{ii},'h'];

        % Construction of names to plot the solution in MATLAB
        Xh{1}.StrPlot{end+1} = SolutionInfo.outputnames{ii};

        % Construction of names to save the CSV file
        if SolutionInfo.outputsizes{ii} == 1
            Xh{1}.StrCSV{end+1}  = [SolutionInfo.outputnames{ii},'h'];
        else
            Xh{1}.StrCSV{end+1}  = [SolutionInfo.outputnames{ii},'h1'];
            Xh{1}.StrCSV{end+1}  = [SolutionInfo.outputnames{ii},'h2'];
        end

        % Construction of the final solution
        Xh{1}.Solution{end+1} = Solh{ii}(1:CountID,:);
    end

    %% Construction of the output exact solution structure
    Xexact = {};

    if Data.PlotExact
        Xexact{1}.StrPlot  = {'x', 'y'};

        Xexact{1}.Solution = { X(1:CountID),   Y(1:CountID)};

        for ii = 1:length(SolutionInfo.outputnames)

            % Construction of names to plot the solution in MATLAB
            Xexact{1}.StrPlot{end+1} = SolutionInfo.outputnames{ii};

            % Construction of the final solution
            Xexact{1}.Solution{end+1} = Solex{ii}(1:CountID,:);
        end
    end


end
