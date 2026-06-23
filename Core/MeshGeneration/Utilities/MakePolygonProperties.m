%> @file  MakePolygonProperties.m
%> @author Mattia Corti
%> @date 22 May 2025
%> @brief Construct area, bounding box and kb associated to each polygonal
%> element
%>
%==========================================================================
%> @section classMakePolygonProperties Class description
%==========================================================================
%> @brief            Construct area, bounding box and kb associated to each polygonal
%> element
%>
%> @param region     Struct containing polygonal mesh
%>
%> @retval region    Struct containing polygonal mesh and the specific
%> polygon properties for each element 
%>
%==========================================================================

function [region] = MakePolygonProperties(region)

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
        region.area(ii,1) = polyarea(region.coords_element{ii}(:,1),region.coords_element{ii}(:,2));

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
end