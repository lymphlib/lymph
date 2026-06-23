%> @file  MetisMextoRegion.m
%> @author Mattia Corti
%> @date 16 May 2023
%> @brief Conversion from metismex mesh format to the lymph one
%>
%==========================================================================
%> @section classMetisMextoRegion Class description
%==========================================================================
%> @brief            Conversion from MetisMex mesh format to the Lymph one.
%> The final mesh is saved in a file with the region and neighbor
%> structures.
%>
%> @param mesh       Mesh in metismex format
%> @param Data       Struct with problem's data
%> @param filename   Name of the file where the final lymph mesh is saved.
%>
%>
%==========================================================================


function MetisMextoRegion(mesh, Data, filename)

    %% Create region.coord

    region.coord = mesh.elem{1};
    
    %% Create region.ne

    region.ne = length(mesh.elem{3});
    
    %% Create region.nedges

    region.nedges = zeros(1,length(mesh.elem{3}));

    for jj = 1:length(mesh.elem{3})
        region.nedges(jj) = length(mesh.elem{3}{jj});
    end
    

    %% Create region.coords_element and region.connectivity

    for kk = 1:length(mesh.elem{3})

        % Extraction of points of the element
        points = mesh.elem{2}(mesh.elem{3}{kk},:);
        points_ord = points(1,:);
        points(1,:) = [];
        
        % Cycle to reorder the points of the mesh element
        while size(points,1)>0
            app2 = points(sum(points == points_ord(end),2)==1,:);
            points_ord = [points_ord, sum(app2.*(1-(app2==points_ord(end))))];
            points(sum(points == points_ord(end-1),2)==1,:) = [];
        end

        % Elimination of the repeated element in the array
        points_ord(end) = [];

        % Extraction of the element coordinates for the kk-th element
        region.coords_element{kk} = mesh.elem{1}(points_ord,:);
        
        % Construction of the region connectivity of the kk-th element
        region.connectivity{kk} = points_ord;
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

    %% Mesh information saving

    save(filename,"region","neighbor")

end