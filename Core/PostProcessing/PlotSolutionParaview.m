%> @file  PlotSolutionParaview.m
%> @author Ilario Mazzieri, Stefano Bonetti, Mattia Corti
%> @date 24 July 2024
%> @brief  Plot of a vector field
%>
%==========================================================================
%> @section classPlotSolutionParaviewLaplacian Class description
%> @brief  Plot of a vector field
%>
%> @param Setup      Struct with problem's setup
%> @param Data       Struct with problem's data
%> @param DispVTK    Cell array with solution to be plotted. It contains all the solutions of the PDE 
%> @param counter    counter id for saving plots
%> @param stringVTK  vector of strings to write in the title (one for each solution of the PDE)
%> @param Id         vector of cell ids associated to the points DispVTK
%> @param Bd         vector of booleans: 1: point on the edge - 0: inner point
%> @param jj         index of the current physics: useful to distinguish variables in multi-physics problems
%>
%> @retval []
%>
%==========================================================================

function [] = PlotSolutionParaview(Setup, Data, DispVTK, counter, stringVTK, Id, Bd, jj)

% Build the delaunay mesh for Paraview
tri = cell(size(Id,1),1);

for ii = 1:size(Id,1)
    Constraint = [Bd{ii} Bd{ii}+1];
    Constraint(end,end) = Constraint(1,1);
    triapp = delaunayTriangulation(DispVTK{1}(Id(ii,1):Id(ii,2)),DispVTK{2}(Id(ii,1):Id(ii,2)),Constraint);
    if ii == 1
        tri{ii} = triapp(isInterior(triapp),:);
    else
        tri{ii} = triapp(isInterior(triapp),:) + max(max(tri{ii-1}));
    end
end
tri = cell2mat(tri);

% Create file name
fname = fullfile(Setup.OutFolderVTK, [Data.name,'_Phys_', num2str(jj), '_', num2str(counter),'.vtk']);

% Write the first solution-field in the .vtk file
WriteVtk(fname, DispVTK{1}, DispVTK{2}, 0*DispVTK{2}, DispVTK{3}, tri, stringVTK{3}, 0);

% Append the other solution-fields (if present) in the .vtk file 
for i = 4:length(DispVTK)
    WriteVtk(fname, DispVTK{1}, DispVTK{2}, 0*DispVTK{2}, DispVTK{i}, tri, stringVTK{i}, 1);
end
