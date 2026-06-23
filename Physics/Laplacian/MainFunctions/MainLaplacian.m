%> @file  MainLaplacian.m
%> @author The Lymph Team
%> @date 5 June 2026
%> @brief Solution of the Poisson problem with PolydG
%>
%==========================================================================
%> @section classMainLaplacian Class description
%==========================================================================
%> @brief            Solution of the Poisson problem with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation's setup
%>
%> @retval Error     Struct for convergence test
%>
%==========================================================================

function [Error] = MainLaplacian(Data,Setup)

fprintf('\nSolution of the Poisson problem \n');
fprintf('... here display some info about the simulation ... \n');
fprintf('\n------------------------------------------------------------------\n')

%% Mesh setup
[mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data);

% SaveName VTK
Data.name = strcat(Data.name,'_',num2str(femregion.nel));

%% Adaptivity flags preallocation
AdaptFlag = 1;
AdaptIts  = 0;

Matrices = [];
F = [];
Solution.Indicator = [];

%% Adaptivity cycle
while AdaptFlag && AdaptIts <= Data.AdaptIts
    %% Matrix Assembly
    fprintf('\nMatrices computation ... \n');
    tic
    [Matrices] = MatrixAssemblyLaplacian(Data, mesh, femregion, Matrices);
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    
    %% Right-hand side assembly
    
    fprintf('\nComputing RHS ... \n');
    tic
    [F] = ForcingTermAssemblyLaplacian(Data, mesh, femregion, F);
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    
    %% Solving the linear system
    fprintf('\nSolving linear system ... ');
    tic
    Solution.U = Matrices.A \ F.F;
    toc
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
    
    %% Compute the indicator for p-adaptivity
    if Data.Adaptivity
        fprintf('\nComputing adaptive indicator ... ');
        tic
        [Solution.Indicator] = ComputeAdaptiveIndicatorLaplacian(Data, mesh, femregion, Solution);
        toc
        fprintf('\nDone\n')
        fprintf('\n------------------------------------------------------------------\n')

        %% Postprocess solution
        if Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution
            PostProcessSolution(Setup, Data, mesh, femregion, AdaptIts, Solution);
        end
        
    end
    %% p-Adaptivity
    if Data.Adaptivity 
        % Compute new femregion
        if AdaptIts < Data.AdaptIts
            [Data, femregion] = ComputeUpdatedDegrees(Data, mesh, femregion, Solution, AdaptIts);

            Solution.Uold = Solution.U;
            if all(femregion.degree == femregion.degree_old)
                AdaptFlag = 0;
            end
        end
        AdaptIts = AdaptIts + 1; 
    else
        AdaptFlag = 0;
    end

end

%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');
    [Error] = ComputeErrorsLaplacian(Data, mesh, femregion, Solution.U);

    Error.nel = femregion.nel;
    Error.h   = Data.h;
    Error.p   = Data.degree;
    Error.NDoF= femregion.ndof;

    Solution.Cells_L2 = Error.Cells_L2;
    Solution.Cells_dG = Error.Cells_dG;
    
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

%% Postprocess solution
if Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution
    PostProcessSolution(Setup, Data, mesh, femregion, AdaptIts, Solution);  
end


fprintf('\nBye \n');


