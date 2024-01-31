%> @file  MakeNeighborDoubleCircle.m
%> @author Ilario Mazzieri
%> @date 16 April 2023
%> @brief Construction of the neighbor structure.
%>
%==========================================================================
%> @section classMakeNeighborDoubleCircle Class description
%==========================================================================
%> @brief  Build neighbor structure for internal and boundary edges
%>
%> @param Data        Struct with problem's data
%> @param region      Struct for mesh elements
%> @param SimType     String simulation type, used for boundary tag       
%>
%> @retval region     Boundary connectivity and boundary tags
%> @retval neighbor   Neigbor struct having fields
%>                     - nedges(i) num of edges for el. i
%>                     - neigh{i}(j) el-id for neigh. of el. i edge j
%>                     - neighedges{i}(j) edge-id for neigh. of el. i edge j 
%>
%==========================================================================
function [region,neighbor]= MakeNeighborDoubleCircle(Data,region,SimType)

%% Make internal neighbor
[neighbor] = MakeNeighborInternal(region);


%% Compute boundary tags and region connectivity

k = 1;
for i = 1 : region.ne
    for j = 1 : neighbor.nedges(i)
        
        id_edge = neighbor.neigh{i}(j);
        
        if (id_edge == -1)
            B_tag = 2;
            if(j < neighbor.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
          
            if(region.coord(edge(1),1)<= 0)
                B_tag = 3;
            end
            
            
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
        end
    end
end

%% Make boundary neighbor
[neighbor] = MakeNeighborBoundary(Data,region,neighbor,SimType);






