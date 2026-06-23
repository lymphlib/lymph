%> @file  ThetaMethodSolver.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Implementation of theta-method for the FHN problem.
%>
%==========================================================================
%> @section classFHNThetaMethodSolver Class description
%==========================================================================
%> @brief            Implementation of theta-method for the FHN problem.
%>
%> @param Setup             Struct with setup information
%> @param Data              Struct with problem's data
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param mesh              Mesh struct (neighbor and region structs)
%> @param Matrices          Matrices associated to the linear terms
%> @param Solution          Initial condition of the problem
%>
%> @retval Solution         Solution struct at the final time
%> @retval femregion        Finite Element struct at the final time
%> @retval IntError	        Estimate of the errors associated with an integral in time
%>
%==========================================================================

function [Solution, femregion, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, Solution)


    %% Right-hand side assembly
    [F_new] = ForcingTermAssemblyFHN(Data, mesh, femregion, Data.t0, Solution.t_old, []);

    %% Second order extrapolation   
    if Data.theta == 0.5
        [Solution.t_oold]   = ComputeModalSolutionFHN(Data, mesh, femregion, Data.t0-Data.dt);
        Solution.t_oold.u_h = Matrices.M_prj\Solution.t_oold.u_h;
        Solution.t_oold.w_h = Matrices.M_prj\Solution.t_oold.w_h;
    end

    %% Initialization of the integral component of the error
    IntError = 0;
    counter  = 1;
    
    %% Time-loop
    for t = Data.t0+Data.dt:Data.dt:Data.T

        if Data.theta == 0.5
            Solution.t_star.u_h = (3*Solution.t_old.u_h - Solution.t_oold.u_h)./2;
            Solution.t_star.w_h = (3*Solution.t_old.w_h - Solution.t_oold.w_h)./2;
        else
            Solution.t_star = Solution.t_old;
        end
       
        fprintf(['\n - Time: ', num2str(t)]);
        
        %% Update of the forcing term
        F_old = F_new;  

        [F_new] = ForcingTermAssemblyFHN(Data, mesh, femregion, t, Solution.t_star, []);
        
        %% Problem resolution with (possibly) nonlinear solver
        Solution.t = PicardIterations(Data, mesh, femregion, Matrices, Solution.t_old, F_old, F_new, t);

        %% Error computation
        if Setup.isError == 1
            [ErrJIT] = ComputeErrorsFHN(Data, mesh, femregion, Solution.t.u_h, t);
            IntError = IntError + Data.dt*ErrJIT.dG.^2;
        end

        %% Compute the indicator for p-adaptivity
        if (mod(counter,Data.AdaptivityStep)==0) && Data.Adaptivity
            [Data, Solution, femregion, Matrices, F_new] = AdaptivityCycleFHN(Data, mesh, femregion, Solution, Matrices, F_old, F_new, t);
        end

        %% Postprocess solution
        if (mod(counter,Data.VisualizationStep)==0) && (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution)
                PostProcessSolution(Setup, Data, mesh, femregion, counter, Solution.t, t);
        end
    
        Solution.t_old = Solution.t;
        %% Time advancement
        if Data.theta == 0.5
            Solution.t_oold = Solution.t_old;
        end

        counter = counter + 1;
        
    end

    disp('------------------------------------------------------')
    disp(['Ending time: ', num2str(t)]);
    disp('------------------------------------------------------')
    
