function [mesh, femregion, h] = MeshFemregionSetup(Setup, Data, TagEl, LabelEl)

    %% Load Region
    fprintf('\nLoading Region ... \n');
    
    mesh = load(Data.meshfile);
    
    % Compute mesh size
    h = max(sqrt((mesh.region.BBox(:,1)-mesh.region.BBox(:,2)).^2  ...
        + (mesh.region.BBox(:,3)-mesh.region.BBox(:,4)).^2));
    
    fprintf(['Number of Polygonal Elements: ', num2str(mesh.region.ne)]);
    fprintf(['\nMesh size: ', num2str(h)]);
    fprintf('\n\n------------------------------------------------------------------\n')

    % checking tags for elements
    for i = 1 : length(mesh.region.id)
        for k = 1 : length(TagEl)
            for j = 1 : length(TagEl{k})
                if mesh.region.id(i) == TagEl{k}(j)
                    mesh.region.tag(i,1) = LabelEl{k};
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
