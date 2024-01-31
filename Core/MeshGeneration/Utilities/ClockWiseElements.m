%> @file  ClockWiseElements.m
%> @author Mattia Corti
%> @date 16 May 2023
%> @brief Control that the vertices of the mesh elements are stored in
%> clockwise order.
%>
%==========================================================================
%> @section classClockWiseElements Class description
%==========================================================================
%> @brief           Control that the vertices of the mesh elements are stored in
%> clockwise order. If not, the order of the vertices is changed to be
%> consistent with the required description.
%>
%> @param region     Mesh region in lymph format.
%>
%> @retval region    Mesh region in lymph format updated.
%>
%>
%==========================================================================

function [region] = ClockWiseElements(region)

    % Cycle over the mesh elements
    for jj = 1:region.ne

        % Extraction of the element barycentres
        xB = sum(region.coords_element{jj}(:,1))/size(region.coords_element{jj},1);
        yB = sum(region.coords_element{jj}(:,2))/size(region.coords_element{jj},1);
        
        % Extraction of the element's distances (barycentre-node)
        V1 = [region.coords_element{jj} zeros(size(region.coords_element{jj},1),1)]-[xB yB 0].*ones(size(region.coords_element{jj},1),3);
        V2 = [[region.coords_element{jj}(2:end,:);region.coords_element{jj}(1,:)] zeros(size(region.coords_element{jj},1),1)]-[region.coords_element{jj}(1:end,:) zeros(size(region.coords_element{jj},1),1)];
        
        % Cross product
        A = cross(V1,V2,2);

        % Clockwise control
        if A(max(abs(A(:,3)))==abs(A(:,3)),3) < 0
            
            % Element rotation
            region.coords_element{jj} = flip(region.coords_element{jj});
            region.connectivity{jj} = flip(region.connectivity{jj});

            % Element visualization
            fprintf(1,'Element %d coordinates are now stored in clockwise order\n', jj);
            plot(xB,yB,'xr')

        else

            % Element visualization
            plot(xB,yB,'xb')
        
        end
        
        % Element visualization
        if mod(jj,100) == 0
            pause(1e-5)
        end

    end
