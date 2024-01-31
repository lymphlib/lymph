function PlotPolyMesh(region)%, neighbour)
%figure;
% axis equal;
% axis off;
figure
hold on;
for i = 1:size(region.coords_element,2)
    XX = region.coords_element{1,i}(:,1);
    YY = region.coords_element{1,i}(:,2);
    C = ones(size(XX));
    plot([XX ; XX(1)],[YY ; YY(1)],'Color','k','LineWidth',1);
    if(size(region.coords_element,2) < 200)
        text(mean(XX),mean(YY),num2str(i),'Color','r');
    end
end
axis equal;
% xlim([Dati.domain(1) Dati.domain(2)]);
% ylim([Dati.domain(3) Dati.domain(4)]);

% for i = 1:size(region.coords_element,2)
%     X = [];
%     Y = [];
%     for j = 1:size(region.connectivity{i})
%
%         %interface edges
%         if neighbour.neigh{i}(j) > 0 && region.id(neighbour.neigh{i}(j)) == 3 && region.id(i) == 1
%             if j < size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','g','LineWidth',1.5)
%             elseif j == size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
%                     [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','g','LineWidth',2)
%             end
%         end
%
%         % Mixed -- one Neumann one Dirichlet
%         if neighbour.neigh{i}(j) == -4
%             if j < size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','m','LineWidth',1.5)
%             elseif j == size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
%                     [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','m','LineWidth',1.5)
%             end
%         end
%         % Absorbing
%         if neighbour.neigh{i}(j) == -3
%             if j < size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','g','LineWidth',1.5)
%             elseif j == size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
%                     [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','g','LineWidth',1.5)
%             end
%         end
%         % Neumann
%         if neighbour.neigh{i}(j) == -2
%             if j < size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','r','LineWidth',1.5)
%             elseif j == size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
%                     [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','r','LineWidth',1.5)
%             end
%         end
%         % Dirichlet
%         if neighbour.neigh{i}(j) == -1
%             if j < size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(j,1),region.coords_element{i}(j+1,1)],[region.coords_element{i}(j,2),region.coords_element{i}(j+1,2)],'Color','b','LineWidth',1.5)
%             elseif j == size(region.connectivity{i},1)
%                 plot([region.coords_element{i}(size(region.connectivity{i},1),1),region.coords_element{i}(1,1)],...
%                     [region.coords_element{i}(size(region.connectivity{i},1),2),region.coords_element{i}(1,2)],'Color','b','LineWidth',1.5)
%             end
%         end
%     end
% end
