%> @file  ThetaMethodSolver.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 15 April 2026
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
%> @retval Solution         Struct containing the problem solution
%> @retval IntError         Struct containing the time integrals of dG and
%> L2 errors
%>
%==========================================================================

function [Solution, femregion, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, u_old)

    Solution.u_old = u_old;
    %% Assembling the constant component of the matrices

    LHS = Matrices.Mprj + Data.dt*Data.theta*(Matrices.A+Matrices.M);
    RHS = Matrices.Mprj - Data.dt*(1-Data.theta)*(Matrices.A+Matrices.M);

    %% Right-hand side assembly
    fprintf('\nComputing forcing terms ... \n');

    [F_new] = ForcingTermAssemblyHeat(Data, mesh, femregion, Data.t0, []);

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

        [F_new] = ForcingTermAssemblyHeat(Data, mesh, femregion, t, []);

        %% Assembling the dynamic component of the matrices in the two treatments of nonlinear term

        % Construction of the complete RHS for the theta method
        F = RHS * Solution.u_old + Data.dt*Data.theta*F_new.F + Data.dt*(1-Data.theta)*F_old.F;

        % Linear system resolution
        Solution.u_h = LHS\F;

        %% Error computation
        if Setup.isError == 1
            [ErrJIT] = ComputeErrorsHeat(Data, mesh, femregion, Solution.u_h, t);
            IntError = IntError + Data.dt*ErrJIT.dG.^2;
        end

        %% Compute the indicator for p-adaptivity
        if (mod(counter,Data.AdaptivityStep)==0) && Data.Adaptivity
            [Data, Solution, femregion, Matrices, RHS, LHS, F_new] = AdaptivityCycle(Data, mesh, femregion, Solution, Matrices, F_old, F_new, t);
        end

        %% Postprocess solution
        if (mod(counter,Data.VisualizationStep)==0) && (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution)
            PostProcessSolution(Setup, Data, mesh, femregion, counter, Solution, t);
        end
        %% Time advancement
        Solution.u_old = Solution.u_h;

        counter = counter + 1;

    end

    disp('------------------------------------------------------')
    disp(['Ending time: ', num2str(t)]);
    disp('------------------------------------------------------')

