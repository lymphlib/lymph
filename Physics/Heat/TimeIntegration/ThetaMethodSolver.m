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

    counter = 1;

    %% Time-loop
    for t = Data.dt:Data.dt:Data.T

        fprintf(['\n - Time: ', num2str(t)]);

        %% Update of the forcing term if needed

        F_old = F_new;

        if ~Data.homog_source_f || Data.TagApplyBCs == 1
            [F_new] = ForcingHeat(Data, mesh.neighbor, femregion, t);
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

        %% Postprocess solution
        if (mod(counter,Data.VisualizationStep)==0) && (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution)
            PostProcessSolution(Setup, Data, mesh, femregion, counter, u_h,t);
        end
        %% Time advancement
        u_old = u_h;

        counter = counter + 1;

    end

    disp('------------------------------------------------------')
    disp(['Ending time: ', num2str(t)]);
    disp('------------------------------------------------------')

