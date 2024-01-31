%> @file  PlotSolutionParaview.m
%> @author Mattia Corti
%> @date 3 October 2023
%> @brief  Plot of a scalar field
%>
%==========================================================================
%> @section classPlotSolutionParaviewHeat Class description
%> @brief  Plot of a scalar field
%
%> @param Setup      Struct with problem's setup
%> @param Data       Struct with problem's data
%> @param Gu         Vector with solution to be plotted
%> @param counter    counter id for saving plots
%>
%> @retval []
%>
%==========================================================================

function PlotSolutionParaview(Setup, Data, Gu, counter)

    if ~exist(Setup.OutFolderVTK,'dir')
        mkdir(Setup.OutFolderVTK); 
    end
    
    % Create the delaunay triangulation
    tri = delaunay(Gu(:,1),Gu(:,2));
    
    fname = fullfile(Setup.OutFolderVTK, [Data.name,'_u_', num2str(counter),'.vtk']);
    WriteVtk(fname, Gu(:,1), Gu(:,2), 0*Gu(:,2), Gu(:,3), tri, 'u');

