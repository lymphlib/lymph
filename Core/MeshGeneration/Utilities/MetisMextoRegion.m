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

        % Element visualization
        plot(mesh.elem{1}(points_ord,1),mesh.elem{1}(points_ord,2))
        hold on

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

    region.max_kb = cell(1,length(region.connectivity));

    % Visualization of computational progress
    prog = 0;
    fprintf(1,'Computation Progress of max kb calculation: %3d%%\n',prog);
    
    % Cycle over the ii-th element
    for ii = 1:length(region.connectivity)
        
        prog = 100*ii/region.ne;
        fprintf(1,'\b\b\b\b%3.0f%%',prog);

        % Extraction of necessary information
        nedge = length(region.connectivity{ii});
        coords_element = region.coord(region.connectivity{ii},:);
        
        % Computation of bounding box
        region.BBox(ii,1) = min(region.coords_element{ii}(:,1));
        region.BBox(ii,2) = max(region.coords_element{ii}(:,1));
        region.BBox(ii,3) = min(region.coords_element{ii}(:,2));
        region.BBox(ii,4) = max(region.coords_element{ii}(:,2));

        % Computation of element area
        region.area(ii) = polyarea(region.coords_element{ii}(:,1),region.coords_element{ii}(:,2));

        % Preallocation of max_kb structure for each element
        region.max_kb{ii} = zeros(nedge,1);
        
        % Cycle over the jj-th point of the element
        for jj = 1:nedge

            % Extraction of the first point
            v(1,:) = coords_element(jj,:);

            % Extraction of the second point
            if jj == nedge
                v(2,:) = coords_element(1,:);
            else
                v(2,:) = coords_element(jj+1,:);
            end

            % Cycle over all the others points of the element
            for kk = 1:nedge

                if (kk ~= jj) && (kk ~= jj+1)
                    
                    % Extraction of the third point
                    v(3,:) = coords_element(kk,:);

                    % Construction of the triangle and area computation
                    [x_tria, y_tria] = poly2cw( v(:,1), v(:,2));
                    area_tria = polyarea(x_tria, y_tria);

                    % Intersect the triangle and the mesh element
                    [x1,y1] = polybool('intersection', coords_element(end:-1:1,1), coords_element(end:-1:1,2), x_tria, y_tria);
                    
                    % Control of the correct construction of the intersection and update of max_kb
                    if (~any(isnan(x1))) && (abs(polyarea(x1,y1) - area_tria) < 1e-6)
                    
                        region.max_kb{ii}(jj) = max(area_tria, region.max_kb{ii}(jj));
                    
                    end
                end
            end
        end
    end

    %% Create region.id

    region.id = ones(size(region.BBox,1),1);
    
    
    %% Create neighbor structure

    [region, neighbor]= MakeNeighbor(Data,region,'laplacian');
    

    %% Create region.tag

    region.tag = 'P';


    %% Mesh information saving

    save(filename,"region","neighbor")

end