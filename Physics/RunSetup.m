%> @file  RunSetup.m
%> @author The Lymph Team
%> @date 9 October 2023
%> @brief Configuration setup
%>
%==========================================================================
%> @section classRunSetup Class description
%==========================================================================
%> @brief            Configuration setup
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================

%% Simulation - Setup

% Plot polygonal mesh y/n
Setup.isPlotMesh = 0;

% Save VTK polygonal mesh y/n
Setup.isSaveVTKMesh = 1;

% Plot solution y/n
Setup.isPlotSolution = 1;

% Save solution y/n -> .mat file
Setup.isSaveSolution = 0;
Setup.OutFolder      = 'Output';

% Additional solution output y/n -> .csv file
Setup.isSaveCSV = 0;

% Additional solution output y/n -> .vtk file
Setup.isSaveVTK = 0;
Setup.OutFolderVTK = 'OutputVTK';

% Compute errors y/n
Setup.isError = 1;
