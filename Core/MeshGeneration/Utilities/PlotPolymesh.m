%> @file  PlotPolymesh.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Plot of the polygonal mesh
%>
%==========================================================================
%> @section classPlotPolymesh Class description
%==========================================================================
%> @brief              Plot of the polygonal mesh and boundary elements
%>
%> @param femregion   Structure containing all the information 
%> about the finite element approximation (See CreateDOF.m)
%> @param neighbor     Structure containing info of neighbor el (See
%MakeNeighbor.m)
%>
%> @retval fig_id      Id of figure
%>
%==========================================================================
function [fig_id] = PlotPolymesh(neighbor, femregion)
    %% Figure initialization
        
    fig_id = figure;
    hold on;


    %% Plot of mesh grid
    
    for i = 1:size(femregion.coords_element,2)
        XX = femregion.coords_element{1,i}(:,1);
        YY = femregion.coords_element{1,i}(:,2);
        plot([XX ; XX(1)],[YY ; YY(1)],'Color','k','LineWidth',1.5);
        if femregion.nel < 100
            text(mean(XX),mean(YY),num2str(i),'Color','r');
        end
    end

    axis equal
        

    %% Plot of boundaries
    
    for i = 1:size(femregion.coords_element,2)
    
        for j = 1:size(femregion.connectivity{i},2)
            
            
            %% Plot of Mixed boundaries 

            if neighbor.neigh{i}(j) == -4
                if j < size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(j,1),femregion.coords_element{i}(j+1,1)],[femregion.coords_element{i}(j,2),femregion.coords_element{i}(j+1,2)],'Color','m','LineWidth',1.5)
                elseif j == size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(size(femregion.connectivity{i},2),1),femregion.coords_element{i}(1,1)],...
                        [femregion.coords_element{i}(size(femregion.connectivity{i},2),2),femregion.coords_element{i}(1,2)],'Color','m','LineWidth',1.5)
                end
            end


            %% Plot of Absorbing boundaries
            
            if neighbor.neigh{i}(j) == -3
                if j < size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(j,1),femregion.coords_element{i}(j+1,1)],[femregion.coords_element{i}(j,2),femregion.coords_element{i}(j+1,2)],'Color','g','LineWidth',1.5)
                elseif j == size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(size(femregion.connectivity{i},2),1),femregion.coords_element{i}(1,1)],...
                        [femregion.coords_element{i}(size(femregion.connectivity{i},2),2),femregion.coords_element{i}(1,2)],'Color','g','LineWidth',1.5)
                end
            end


            %% Plot of Neumann boundaries
            
            if neighbor.neigh{i}(j) == -2
                if j < size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(j,1),femregion.coords_element{i}(j+1,1)],[femregion.coords_element{i}(j,2),femregion.coords_element{i}(j+1,2)],'Color','r','LineWidth',1.5)
                elseif j == size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(size(femregion.connectivity{i},2),1),femregion.coords_element{i}(1,1)],...
                        [femregion.coords_element{i}(size(femregion.connectivity{i},2),2),femregion.coords_element{i}(1,2)],'Color','r','LineWidth',1.5)
                end
            end


            %% Plot of Dirichlet boundaries
            
            if neighbor.neigh{i}(j) == -1
                if j < size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(j,1),femregion.coords_element{i}(j+1,1)],[femregion.coords_element{i}(j,2),femregion.coords_element{i}(j+1,2)],'Color','b','LineWidth',1.5)
                elseif j == size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(size(femregion.connectivity{i},2),1),femregion.coords_element{i}(1,1)],...
                        [femregion.coords_element{i}(size(femregion.connectivity{i},2),2),femregion.coords_element{i}(1,2)],'Color','b','LineWidth',1.5)
                end
            end

            %% Plot of interface boundaries

            if neighbor.neigh{i}(j)>0 && femregion.tag(neighbor.neigh{i}(j)) ~= femregion.tag(i)
                if j < size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(j,1),femregion.coords_element{i}(j+1,1)],[femregion.coords_element{i}(j,2),femregion.coords_element{i}(j+1,2)],'Color','y','LineWidth',1.5)
                elseif j == size(femregion.connectivity{i},2)
                    plot([femregion.coords_element{i}(size(femregion.connectivity{i},2),1),femregion.coords_element{i}(1,1)],...
                        [femregion.coords_element{i}(size(femregion.connectivity{i},2),2),femregion.coords_element{i}(1,2)],'Color','y','LineWidth',2)
                end
            end

        end

    end
    
end
