%> @file  ThetaMethodSolver.m
%> @author Mattia Corti
%> @date 11 November 2023
%> @brief Implementation of theta-method for the Fisher-KPP problem.
%>
%==========================================================================
%> @section classFisherKPPThetaMethodSolver Class description
%==========================================================================
%> @brief            Implementation of theta-method for the Fisher-KPP problem.
%>
%> @param Setup             Struct with setup information
%> @param Data              Struct with problem's data
%> @param femregion         Finite Element struct (see CreateDOF.m)
%> @param mesh              Mesh struct (neighbor and region structs)
%> @param Matrices          Matrices associated to the linear terms
%> @param c_old             Initial condition of the problem (t=0)
%>
%> @retval Solution         Solution at the final time
%> @retval femregion        Finite Element struct at the final time
%> @retval IntError	    Estimate of the errors associated with an integral in time
%>
%==========================================================================

function [Solution, femregion, IntError] = ThetaMethodSolver(Setup, Data, femregion, mesh, Matrices, c_old)

    Solution.c_old = c_old;
    
    %% Right-hand side assembly
    [F_new] = ForcingTermAssemblyFKPP(Data, mesh, femregion, Data.t0, []);

    %% Second order extrapolation   
    if Data.theta == 0.5
        Solution.c_oold = ComputeModalSolutionFKPP(Data, mesh, femregion, Data.t0-Data.dt);
        Solution.c_oold = Matrices.M_prj\Solution.c_oold;
    else
        Solution.c_oold = Solution.c_old;
    end

    %% Initialization of the integral component of the error
    IntError = 0;
    counter  = 1;
    %% Time-loop
    for t = Data.t0+Data.dt:Data.dt:Data.T
       
        fprintf(['\n - Time: ', num2str(t)]);
        
        %% Update of the forcing term
        F_old = F_new;  

        [F_new] = ForcingTermAssemblyFKPP(Data, mesh, femregion, t, []);
        
        %% Problem resolution with (possibly) nonlinear solver
        Solution.c_h = PicardIterations(Data, femregion, Matrices, Solution, F_old, F_new);

        %% Error computation
        if Setup.isError == 1
            [ErrJIT] = ComputeErrorsFKPP(Data, mesh, femregion, Solution.c_h, t);
            IntError = IntError + Data.dt*ErrJIT.dG.^2;
        end

        %% Compute the indicator for p-adaptivity
        if (mod(counter,Data.AdaptivityStep)==0) && Data.Adaptivity
            [Data, Solution, femregion, Matrices, F_new] = AdaptivityCycleFKPP(Data, mesh, femregion, Solution, Matrices, F_old, F_new, t);
        end

        %% Postprocess solution
        if (mod(counter,Data.VisualizationStep)==0) && (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution)
                PostProcessSolution(Setup, Data, mesh, femregion, counter, Solution, t);
        end
    
        %% Time advancement
        if Data.theta == 0.5
            Solution.c_oold = Solution.c_old;
        end
        Solution.c_old  = Solution.c_h;

        counter = counter + 1;
        
    end

    disp('------------------------------------------------------')
    disp(['Ending time: ', num2str(t)]);
    disp('------------------------------------------------------')
    
