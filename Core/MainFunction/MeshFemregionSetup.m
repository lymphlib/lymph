%> @file  MeshFemregionSetup.m
%> @author Mattia Corti
%> @date 29 May 2026
%> @brief Load of the mesh and creation of femregion structure
%>
%==========================================================================
%> @section classMeshFemregionSetup Class description
%==========================================================================
%> @brief            Load of the mesh and creation of femregion structure
%
%> @param Setup      Simulation's setup
%> @param Data       Struct with problem's data
%>
%> @retval mesh      Struct containing mesh information
%> @retval femregion Struct containing finite element space information
%> @retval h         Mesh size
%>
%==========================================================================

function [mesh, femregion, h] = MeshFemregionSetup(Setup, Data)

    %% Load Region
    fprintf('\nLoading Region ... \n');
    
    if contains(Data.meshfile,".mat")
        mesh = load(Data.meshfile);
    elseif contains(Data.meshfile,'.vtk')
        mesh = VTKtoRegion(Data);
    end

    % Compute mesh size
    h = max(sqrt((mesh.region.BBox(:,1)-mesh.region.BBox(:,2)).^2  ...
        + (mesh.region.BBox(:,3)-mesh.region.BBox(:,4)).^2));
    
    fprintf(['Number of Polygonal Elements: ', num2str(mesh.region.ne)]);
    fprintf(['\nMesh size: ', num2str(h)]);
    fprintf('\n\n------------------------------------------------------------------\n')

    % checking tags for elements
    for i = 1 : length(mesh.region.id)
        for k = 1 : length(Data.TagEl)
            for j = 1 : length(Data.TagEl{k})
                if mesh.region.id(i) == Data.TagEl{k}(j)
                    mesh.region.label(i,1) = Data.LabEl{k};
                end
            end
        end
    end

    %% Creation of the finite element space
    fprintf('\nMake femregion ... ')
    
    [femregion] = CreateDOF(Data, mesh.region);
    
    fprintf('\nDone\n')
    fprintf('\n------------------------------------------------------------------\n')
       
    %% Plot polygonal mesh
    if Setup.isPlotMesh
        [~] = PlotPolymesh(mesh.neighbor, femregion);
    end
    
    %% Save VTK polygonal mesh
    if Setup.isSaveVTKMesh
        CreatePolygonalVTK(Data, Setup, mesh.region);
    end

end
