%> @file  PlotSolutionParaview.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Plot of a vector field
%>
%==========================================================================
%> @section classPlotSolutionParaviewLaplacian Class description
%> @brief  Plot of a vector field
%>
%> @param Setup      Struct with problem's setup
%> @param Data       Struct with problem's data
%> @param Disp       Matrix with solution to be plotted
%> @param counter    counter id for saving plots
%> @param string     string to write in the title
%>
%> @retval []
%>
%==========================================================================

function [] = PlotSolutionParaview(Setup, Data, Disp, counter, string)


% Elastic displacement
tri = delaunay(Disp(:,1),Disp(:,2));

fname = fullfile(Setup.OutFolderVTK, [Data.name,'_',string,'_', num2str(counter),'.vtk']);
WriteVtk(fname, Disp(:,1), Disp(:,2), 0*Disp(:,2), [Disp(:,3) Disp(:,4)], tri, string);

