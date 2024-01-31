%> @file  MakeNeighborInternal.m
%> @author Ilario Mazzieri & Paola Antonietti
%> @date 16 April 2023
%> @brief Construction of the neighbor structure for internal edges.
%>
%==========================================================================
%> @section classMakeNeighborInternal Class description
%==========================================================================
%> @brief  Build neighbor structure for internal edges
%>
%> @param region      Struct for mesh elements
%>
%> @retval neighbor   Neigbor struct having fields
%>                     - nedges(i) num of edges for el. i
%>                     - neigh{i}(j) el-id for neigh. of el. i edge j
%>                     - neighedges{i}(j) edge-id for neigh. of el. i edge j 
%>
%==========================================================================
function [neighbor] = MakeNeighborInternal(region)

ne = region.ne;
connectivity = region.connectivity;
neigh = cell(1,ne);
neighedges = cell(1,ne);

for i=1:ne
    neigh{i}=-ones(size(connectivity{i}));
    neighedges{i}=-ones(size(connectivity{i}));
end

for i=1:(ne-1)
    edges =[];
    n_edges = length(connectivity{i});
    
    for vertices = 1:n_edges
        v(vertices)=connectivity{i}(vertices);
    end
    
    for e = 1:n_edges-1
        edges(e,:)=[v(e) v(e+1)];
    end
    edges(n_edges,:) = [v(n_edges) v(1)];
    
    for j=(i+1):ne
        edgesn =[];
        n_edgesn = length(connectivity{j});
        for verticesn = 1:n_edgesn
            vn(verticesn)=connectivity{j}(verticesn);
        end
        
        for e = 1:n_edgesn-1
            edgesn(e,:)=[vn(e) vn(e+1)];
        end
        edgesn(n_edgesn,:) = [vn(n_edgesn) vn(1)];
        
        for s = 1:size(edges,1)
            for t = 1:size(edgesn,1)
                if (edges(s,1) == edgesn(t,2) && edges(s,2) == edgesn(t,1))
                    neigh{i}(s)=j;
                    neigh{j}(t)=i;
                    neighedges{i}(s)=t;
                    neighedges{j}(t)=s;
                end
            end
        end
        
    end
    
    
end

%% Output neighbor structure
neighbor.nedges     = region.nedges;
neighbor.neigh      = neigh;
neighbor.neighedges = neighedges;