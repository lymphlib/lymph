%> @file  VTKtoRegion.m
%> @author Mattia Corti
%> @date 22 May 2025
%> @brief Conversion from VTK mesh format to the lymph one
%>
%==========================================================================
%> @section classVTKtoRegion Class description
%==========================================================================
%> @brief            Conversion from VTK mesh format to the Lymph one.
%> The final mesh is saved in a file with the region and neighbor
%> structures.
%>
%> @param Data       Struct with problem's data
%>
%> @retval mesh      Struct containing mesh region and neighbor
%>
%==========================================================================

function [mesh] = VTKtoRegion(Data)

    disp("Conversion from VTK to lymph mesh")
    %% File opening
    file = fopen(Data.meshfile);
    
    %% Find the points of the mesh
    npoints = 0;
    ii = -1;
    while ii < npoints
    
        %% Find the points positioning in the file
        line = fgets(file,1000);
        if contains(line,"POINTS ")
            npoints = erase(line,"POINTS ");
            npoints = erase(npoints," float");
            npoints = str2num(npoints);
            coords = zeros(3*npoints,1);
            ii = 0;
        end
    
        %% Extract the point coordinates
        if npoints > 0 && ~contains(line,"POINTS ")
           points = split(line);
            for jj = 1:size(points,1)-1
                coords(ii*3+jj,:) = str2num(points{jj});
            end
            ii = ii + 3;
        end
    end
    
    %% Reshape and construct the informations for the region struct
    coords = reshape(coords, [3 npoints])';
    coords(:,3) = [];
    
    region.coord = coords;
    
    %% Find the cells of the mesh
    ncells = 0;
    maxpts = 0;
    ii = -1;
    while ii < ncells
    
        %% Find the cells positioning in the file
        line = fgets(file,1000);
        if contains(line,"CELLS ") || contains(line,"POLYGONS ")
            if contains(line,"CELLS ")
                ncells = erase(line,"CELLS ");
            elseif contains(line,"POLYGONS ")
                ncells = erase(line,"POLYGONS ");
            end
            ncellsapp = split(ncells);
            ncells = str2num(ncellsapp{1});
            maxpts = str2num(ncellsapp{2});
            ii = 0;
            line = fgets(file,1000);
            line = fgets(file,1000);
            maxptspoly = zeros(ncells,1);
        end
    
        % Extract the number of points of each polygon
        if maxpts > 0 && ~contains(line,"CELLS ")
           points = split(line);
            for jj = 1:size(points,1)-1
                maxptspoly(ii+1) = str2num(points{jj});
                ii = ii + 1;
            end
        end
    end
    
    region.ne = ncells-1;
    
    %% Compute number of edges for each polygon
    region.nedges = (maxptspoly(2:end) - maxptspoly(1:end-1))';
    
    line = fgets(file,1000);
    ii = 0;
    connectivity = zeros(maxptspoly(end),1);
    
    while ii < maxptspoly(end)
    
        line = fgets(file,1000);
    
        % Extract the connectivity points
        conn = split(line);
        for jj = 1:size(conn,1)-1
            connectivity(ii+1) = str2num(conn{jj})+1;
            ii = ii + 1;
        end
        
    end
    
    %% Construct connectivity and coordinates of single polygon
    
    for ii = 1:length(maxptspoly)-1
        region.connectivity{ii} = connectivity(maxptspoly(ii)+1:maxptspoly(ii+1))';
        region.coords_element{ii} = coords(region.connectivity{ii},:);
    end
    
    %% Clockwise ordering control

    [region] = ClockWiseElements(region);    
    
    %% Create region.max_kb, region.area and region.BBox

    [region] = MakePolygonProperties(region);

    %% Create region.id
    
    region.id = ones(size(region.BBox,1),1);
    
    %% Create neighbor structure
    
    [region, neighbor]= MakeNeighbor(Data,region,'laplacian');
    
    %% Create region.label
    
    region.label = 'P';
    
    %%  Create output struct
    mesh.region = region;
    mesh.neighbor = neighbor;
    
    %% Mesh information saving
    save([Data.meshfile(1:end-4),'.mat'],"region","neighbor")
