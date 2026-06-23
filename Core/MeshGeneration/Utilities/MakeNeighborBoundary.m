%> @file  MakeNeighborBoundary.m
%> @author Ilario Mazzieri
%> @date 29 May 2026
%> @brief Construction of the neighbor structure for boundary edges.
%>
%==========================================================================
%> @section classMakeNeighborBoundary Class description
%==========================================================================
%> @brief  Build neighbor structure for internal edges. The following
%conditions are implemented: [D N A] = [-1 -2 -3]. D = Dirichlet, N =
%Neumann, A = Absorbing.
%>
%> @param Data        Struct with problem's data
%> @param region      Struct for mesh elements
%> @param neighbor    Struct with nighbor data
%>
%> @retval neighbor   Neigbor struct having fields
%>                     - nedges(i) num of edges for el. i
%>                     - neigh{i}(j) el-id for neigh. of el. i edge j
%>                     - neighedges{i}(j) edge-id for neigh. of el. i edge j
%>
%==========================================================================
function [neighbor] = MakeNeighborBoundary(Data, region, neighbor)

for ii = 1:length(Data.TagBc)
    TagBoundary(Data.TagBc{ii}) = Data.LabBc{ii};
end

counter = 0;
for i = 1 : region.ne
    
    for j = 1 : size(neighbor.neigh{i},2)
        
        if (neighbor.neigh{i}(j) < 0)
            counter = counter +1;
            n_edges = length(region.connectivity{i});
            
            for vertices = 1:n_edges
                v(vertices)=region.connectivity{i}(vertices);
            end
            
            for e = 1:n_edges-1
                edges(e,:)=[v(e) v(e+1)];
            end
            edges(n_edges,:) = [v(n_edges) v(1)];
            
            for k = 1 : size(region.connectivity_bc,1)
                if ((region.connectivity_bc(k,1) == edges(j,1) && (region.connectivity_bc(k,2) == edges(j,2))) || ...
                        (region.connectivity_bc(k,2) == edges(j,1) && (region.connectivity_bc(k,1) == edges(j,2))) )
                    
                    switch TagBoundary(region.bc_tag(k))
                        case('D')
                            neighbor.neigh{i}(j) = -1;
                        case('N')
                            neighbor.neigh{i}(j) = -2;
                        case('A')
                            neighbor.neigh{i}(j) = -3;
                        otherwise
                            disp('Bc not known!')
                    end
                end
            end
        end
    end
end
