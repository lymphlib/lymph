%> @file  MakeNeighborQuadridomain.m
%> @author Ilario Mazzieri
%> @date 16 April 2023
%> @brief Construction of the neighbor structure.
%>
%==========================================================================
%> @section classMakeNeighborQuadridomain Class description
%==========================================================================
%> @brief  Build neighbor structure for internal and boundary edges
%
%> @param Data        Struct with problem's data
%> @param region      Struct for mesh elements
%> @param SimType     String simulation type, used for boundary tag       
%
%> @retval region     Boundary connectivity and boundary tags
%> @retval neighbor   Neigbor struct having fields
%>                     - nedges(i) num of edges for el. i
%>                     - neigh{i}(j) el-id for neigh. of el. i edge j
%>                     - neighedges{i}(j) edge-id for neigh. of el. i edge j 
%>
%==========================================================================
function [region,neighbor]= MakeNeighborQuadridomain(Data,region,SimType)

%% Make internal neighbor
[neighbor] = MakeNeighborInternal(region);


%% Compute boundary tags and region connectivity
k = 1;
for i = 1 : region.ne

    for j = 1 : neighbor.nedges(i)
        id_edge = neighbor.neigh{i}(j);
        if (region.id(i) == 1 && id_edge == -1)
            B_tag = 5;
            if(j < neighbor.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
            
            % left boundary
            if(abs(region.coord(edge(1),1) - Data.domain(1)) < 1.e-8 && ...
                    abs(region.coord(edge(2),1) - Data.domain(1)) < 1.e-8)
                B_tag = 12;
            end
            
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
        elseif (region.id(i) == 2 && id_edge == -1)
            B_tag = 11;
            if(j < neighbor.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
            
            % top boundary
            if(abs(region.coord(edge(1),2) - Data.domain(4)) < 1.e-8 && ...
                    abs(region.coord(edge(2),2) - Data.domain(4)) < 1.e-8)
                B_tag = 10;
            end
            
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
        elseif (region.id(i) == 3 && id_edge == -1)
            B_tag = 8;
            if(j < neighbor.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
            
            % top boundary
            if(abs(region.coord(edge(1),2) - Data.domain(4)) < 1.e-8 && ...
                    abs(region.coord(edge(2),2) - Data.domain(4)) < 1.e-8)
                B_tag = 9;
            end
            
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
        elseif (region.id(i) == 4 && id_edge == -1)
            B_tag = 7;
            if(j < neighbor.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
            
            % bottom boundary
            if(abs(region.coord(edge(1),2) - Data.domain(3)) < 1.e-8 && ...
                    abs(region.coord(edge(2),2) - Data.domain(3)) < 1.e-8)
                B_tag = 6;
            end
            
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
        end
        
        
    end
    
end

%% Make boundary neighbor
[neighbor] = MakeNeighborBoundary(Data,region,neighbor,SimType);

