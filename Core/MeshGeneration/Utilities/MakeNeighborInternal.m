%> @file  MakeNeighborInternal.m
%> @author Ilario Mazzieri, Paola Antonietti, Mattia Corti
%> @date 4 October 2024
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

disp("Neighbor structure construction:");

% Extract number of elements and mesh connectivity
ne           = region.ne;
connectivity = region.connectivity;

% Create neighbor structures
neigh      = cell(1,ne);
neighedges = cell(1,ne);
edges      = cell(1,ne);

for ii = 1:ne
    neigh{ii}      = -ones(size(connectivity{ii}));
    neighedges{ii} = -ones(size(connectivity{ii}));

    % Create edges structure
    v         = connectivity{ii};
    edges{ii} = [v; v(2:end), v(1)]';

end

for ii = 1:(ne-1)

    if any(neigh{ii}<0)
        
        edges_ii = edges{ii};

        for jj = (ii+1):ne

            edges_jj = edges{jj};

            if any(connectivity{ii} == connectivity{jj}','all')
                for ss = 1:size(edges_ii,1)
                    % Control only if not already assigned
                    for tt = 1:size(edges_jj,1)
                        if edges_ii(ss,1) == edges_jj(tt,2)
                            if edges_ii(ss,2) == edges_jj(tt,1)
                                neigh{ii}(ss) = jj;
                                neigh{jj}(tt) = ii;
                                neighedges{ii}(ss) = tt;
                                neighedges{jj}(tt) = ss;
                            end
                        end
                    end
                end
            end
        end

    end


    if mod(ii,1000) == 0
        strprint = strcat(" - ", strcat(num2str(ii))," elements processed!");
        disp(strprint)
    end
    
end

%% Output neighbor structure
neighbor.nedges     = region.nedges;
neighbor.neigh      = neigh;
neighbor.neighedges = neighedges;
