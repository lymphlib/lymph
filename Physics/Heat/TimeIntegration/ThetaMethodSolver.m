%> @file  ThetaMethodSolver.m
%> @author Mattia Corti
%> @date 3 October 2023
%> @brief Implementation of theta-method for the heat equation.
%>
%==========================================================================
%> @section classThetaMethodSolver Class description
%==========================================================================
%> @brief            Implementation of theta-method for the heat equation.
%
%> @param Setup             Struct with setup information
%> @param Data              Struct with problem's data
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param mesh              Mesh struct (neighbor and region structs)
%> @param Matrices          Matrices associated to the PDE
%> @param u_old             Initial condition of the problem (t=0)
%>
%==========================================================================

function [u_h, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, u_old)


    %% Assembling the constant component of the matrices

    LHS = Matrices.M_prj + Data.dt*Data.theta*(Matrices.A+Matrices.M);
    RHS = Matrices.M_prj - Data.dt*(1-Data.theta)*(Matrices.A+Matrices.M);

    %% Right-hand side assembly

    fprintf('\nComputing forcing terms ... \n');

    [F_new] = ForcingHeat(Data, mesh.neighbor, femregion, Data.t0);

    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')

    % Initialization of the integral component of the error
    IntError = 0;

    %% Time-loop
    for t = Data.dt:Data.dt:Data.T

        fprintf(['\n - Time: ', num2str(t)]);

        %% Update of the forcing term if needed

        F_old = F_new;

        if ~Data.homog_source_f || Data.TagApplyBCs == 1

            fprintf('\nComputing forcing terms ... \n');

            [F_new] = ForcingHeat(Data, mesh.neighbor, femregion, t);

            fprintf('\nDone\n')
            fprintf('\n------------------------------------------------------------------\n')

        end


        %% Assembling the dynamic component of the matrices in the two treatments of nonlinear term

        % Construction of the complete RHS for the theta method
        F = RHS * u_old + Data.dt*Data.theta*F_new + Data.dt*(1-Data.theta)*F_old;

        % Linear system resolution
        u_h = LHS\F;

        if Setup.isError == 1
            [ErrJIT] = ComputeErrors(Data, Matrices, femregion, u_h, t);
            IntError = IntError + Data.dt*ErrJIT.err_u_dG.^2;
        end

        %% Time advancement
        u_old = u_h;

        %% Plot of the solution      
        if Setup.isPlotSolution == 1 && mod(t,Data.VisualizationStep) < 0.01*Data.dt
            [Gu] = GetSolutionQuadPoints(Data, femregion, mesh.neighbor, u_old, t);
            ScatterPlot(Gu, mesh.region, Data, 'u');
        end

        %% Save the numerical solution
        if Setup.isSaveSolution == 1 && mod(t, Data.SaveSolutionStep) < 0.01*Data.dt
            filenameout = fullfile(Setup.OutFolder,['Solution_', num2str(ceil(t/Data.dt)), '.mat']);
            save(filenameout);
        end

        %% Save solutions paraview
        if Setup.isSaveCSV && mod(t, Data.SaveSolutionStep) < 0.01*Data.dt
            [Gu] = GetSolutionQuadPoints(Data, femregion, mesh.neighbor, u_old, t);
            filenameout = fullfile(Setup.OutFolder,['Solution_', num2str(ceil(t/Data.dt)), '.csv']);
            
            if Data.PlotExact
                title = {'x','y','u_h','u_ex'};
            else
                title = {'x','y','u_h'};
            end
            Gu = array2table(Gu);
            Gu.Properties.VariableNames(1:length(title)) = title;
            writetable(Gu,filenameout);
        end
        
        if Setup.isSaveVTK && mod(t, Data.SaveSolutionStep) < 0.01*Data.dt
            [Gu] = GetSolutionQuadPoints(Data, femregion, mesh.neighbor, u_old, t);
            PlotSolutionParaview(Setup, Data, Gu, round(t/Data.dt));
        end

    end

    disp('------------------------------------------------------')
    disp(['Ending time: ', num2str(t)]);
    disp('------------------------------------------------------')

