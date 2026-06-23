%> @file  PostProcessSolution.m
%> @author Mattia Corti, Stefano Bonetti
%> @date 11 May 2026
%> @brief Postprocess the numerical solution
%>
%==========================================================================
%> @section classPostProcessSolution Class description
%==========================================================================
%> @brief            Postprocess the numerical solution
%>
%> @param Setup      Struct with setup properties of the problem
%> @param Data       Struct with problem's data
%> @param mesh       Struct with the mesh information
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param counter    Auxiliary index for saving
%> @param Solution   Numerical solution
%> @param time       Current time
%> @param par        Parameter used to choose the background call of the
%> postprocessing procedure
%>
%==========================================================================

function PostProcessSolution(Setup, Data, mesh, femregion, counter, Solution, time, par)

    if ~Setup.isParallel || nargin == 8
        %% Creation of output folders
        if ~exist(Setup.OutFolder,'dir')
            mkdir(Setup.OutFolder);
        end
        
        if ~exist(Setup.OutFolderVTK,'dir')
            mkdir(Setup.OutFolderVTK)
        end
        
        %% Save .mat file
        if Setup.isSaveSolution
            filename = fullfile(Setup.OutFolder,[Data.name,'_',num2str(counter),'.mat']);
            save(filename);
        end
        
        %% Save solutions paraview
        if (Setup.isSaveCSV || Setup.isSaveVTK || Setup.isPlotSolution)
            fprintf('\nGetting solution for post-processing... \n')
            tic
            if nargin == 6
                [XhPlot, XexactPlot] = GetSolutionPostProcessing(Data, femregion, mesh.neighbor, Solution);
            else
                [XhPlot, XexactPlot] = GetSolutionPostProcessing(Data, femregion, mesh.neighbor, Solution, time);
            end
            toc
            fprintf('\nDone\n')
            fprintf('\n------------------------------------------------------------------\n')
        end
        
        %% Save .csv file
        if Setup.isSaveCSV
            fprintf('\nSaving CSV ... ')
            tic
            for jj = 1:length(XhPlot)
                filename = ['Solution_Phys_', num2str(jj),'_',num2str(counter),'.csv'];
                filenameout = fullfile(Setup.OutFolder, filename);
                XhPlotTab = array2table(cell2mat(XhPlot{jj}.Solution), 'VariableNames',XhPlot{jj}.StrCSV);
                writetable(XhPlotTab,filenameout);
            end
            toc
            fprintf('\nDone\n')
            fprintf('\n------------------------------------------------------------------\n')
        end
        
        %% Save .vtk file
        if Setup.isSaveVTK
            fprintf('\nSaving VTK ... ')
            tic
            for jj = 1:length(XhPlot)
                PlotSolutionParaview(Setup, Data, XhPlot{jj}.Solution, counter, XhPlot{jj}.StrVTK, XhPlot{jj}.Id, XhPlot{jj}.Bd, jj);
            end
            toc
            fprintf('\nDone\n')
            fprintf('\n------------------------------------------------------------------\n')
        end
        
        %% Plot solution matlab
        if Setup.isPlotSolution
            close all
            for jj = 1:length(XhPlot)
                if Data.PlotExact
                    ScatterPlot(XhPlot{jj}, XexactPlot{jj}, mesh.region, Data);
                else
                    ScatterPlot(XhPlot{jj}, [], mesh.region, Data);
                end
                
                if isfield(Data,"Adaptivity")
                    if Data.Adaptivity
                        ScatterPlotIndicator(XhPlot{jj}, mesh.region, Data);
                    end
                end
            end
        
        end
    else
        %% Run the postprocessing in background
        parfeval(@PostProcessSolution, 0, Setup, Data, mesh, femregion, counter, Solution, time, 1);
    end
end
