%> @file  ScatterPlot.m
%> @author Stefano Bonetti, Mattia Corti, Ilario Mazzieri
%> @date 27 July 2024
%> @brief Scatter visualization of the elastic solution
%>
%==========================================================================
%> @section classScatterPlot Class description
%==========================================================================
%> @brief            Scatter visualization of the elastic solution
%
%> @param Xh        Cell array that contains the approximate solutions. It has the following form:
%> Xh = { x , y , u_h (x,y), v_h(x,y), ... } for scalar-valued outputs
%> Xh = { x , y , [(u_h)_x(x,y) (u_h)_y(x,y)], [(v_h)_x(x,y) (v_h)_y(x,y)],
%... } for vector-valued outputs
%> @param Xexact    Cell array that contains the exact solutions. It has the following form:
%> Xexact = { x , y , u_ex (x,y), v_ex(x,y), ... } for scalar-valued outputs
%> Xexact = { x , y , [(u_ex)_x(x,y) (u_ex)_y(x,y)], [(v_ex)_x(x,y) (v_ex)_y(x,y)],
%... } for vector-valued outputs
%> where x, y are the coordinates of the quadrature nodes,
%> @param region     Mesh region struct
%> @param Data       Struct with problem's data
%
%> @retval ~
%>
%==========================================================================

function ScatterPlot(Xh, Xex, region, Data)

for i = 3:length(Xh.Solution)

    %% Check if data is scalar or vectorial
    if Data.PlotExact
        if size(Xh.Solution{i},2)~=size(Xex.Solution{i},2)
            error('ScatterPlot: approximate and exact solutions must have the same dimensions');
        end
    end

    if ~(size(Xh.Solution{i},2)==1 || size(Xh.Solution{i},2)==2)
        error(strcat('ScatterPlot: arguments to plot must have either 1 or 2 columns (found ', size(Xh{i},2),')'));
    end

    if size(Xh.Solution{i},2) == 1 % plot of a scalar-variable
    
        figure
 
        %% Control if plot the exact solution
        if Data.PlotExact

            %%  dG-solution scatter plot
            subplot(1,3,1)
            scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}, 'filled');
            title(strcat({'$'},{Xh.StrPlot{i}},{'^\mathrm{h}$'}),'Interpreter','latex','Fontsize',14);
            colorbar
            hold on
            if Data.PlotGridSol
                for kk = 1:region.ne
                    plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                end
            end
            axis equal
            xlim([Data.domain(1), Data.domain(2)])
            ylim([Data.domain(3), Data.domain(4)])

            %% Exact-solution scatter plot
            subplot(1,3,2)
            scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xex.Solution{i}, 'filled');
            title(strcat({'$'},{Xh.StrPlot{i}},{'^\mathrm{ex}$'}),'Interpreter','latex','Fontsize',14);
            colorbar
            hold on
            if Data.PlotGridSol
                for kk = 1:region.ne
                    plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                end
            end
            axis equal
            xlim([Data.domain(1), Data.domain(2)])
            ylim([Data.domain(3), Data.domain(4)])
            
            %% Error between exact and numerical solution scatter plot
            subplot(1,3,3)
            scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}-Xex.Solution{i}, 'filled');
            title(strcat({'$'},{Xh.StrPlot{i}},{'^\mathrm{h}-'},{Xh.StrPlot{i}},{'^\mathrm{ex}$'}),'Interpreter','latex','Fontsize',14);
            colorbar
            hold on
            if Data.PlotGridSol
                for kk = 1:region.ne
                    plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                end
            end
            axis equal
            xlim([Data.domain(1), Data.domain(2)])
            ylim([Data.domain(3), Data.domain(4)])
            
        else

            %%  dG-solution scatter plot
            scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}, 'filled');
            title(strcat({'$'},{Xh.StrPlot{i}},{'^\mathrm{h}$'}),'Interpreter','latex','Fontsize',14);
            colorbar
            hold on
            if Data.PlotGridSol
                for kk = 1:region.ne
                    plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                end
            end
            axis equal
            xlim([Data.domain(1), Data.domain(2)])
            ylim([Data.domain(3), Data.domain(4)])

        end
                
    else % plot of a vector-variable

        for Comp = 1:2 % loop on the two components of the vector-field
            
            figure
        
            %% Control if plot the exact solution
            if Data.PlotExact

                %%  DG-solution scatter plot
                subplot(1,3,1)
                if Comp == 1
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}(:,1), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{h})_x$'}),'Interpreter','latex','Fontsize',14);
                else
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}(:,2), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{h})_y$'}),'Interpreter','latex','Fontsize',14);
                end
                colorbar
                hold on
                if Data.PlotGridSol
                    for kk = 1:region.ne
                        plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                    end
                end
                axis equal
                xlim([Data.domain(1), Data.domain(2)])
                ylim([Data.domain(3), Data.domain(4)])

                %% Exact-solution scatter plot
                subplot(1,3,2)
                if Comp == 1
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xex.Solution{i}(:,1), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{ex})_x$'}),'Interpreter','latex','Fontsize',14);
                else
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xex.Solution{i}(:,2), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{ex})_y$'}),'Interpreter','latex','Fontsize',14);
                end
                colorbar
                hold on
                if Data.PlotGridSol
                    for kk = 1:region.ne
                        plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                    end
                end
                axis equal
                xlim([Data.domain(1), Data.domain(2)])
                ylim([Data.domain(3), Data.domain(4)])
                
                %% Error between exact and numerical solution scatter plot
                subplot(1,3,3)
                if Comp == 1
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}(:,1) - Xex.Solution{i}(:,1), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{h})_x-('},{Xh.StrPlot{i}},{'^\mathrm{ex})_x$'}),'Interpreter','latex','Fontsize',14);
                else
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}(:,2) - Xex.Solution{i}(:,2), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{h})_y-('},{Xh.StrPlot{i}},{'^\mathrm{ex})_y$'}),'Interpreter','latex','Fontsize',14);
                end
                colorbar
                hold on
                if Data.PlotGridSol
                    for kk = 1:region.ne
                        plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                    end
                end
                axis equal
                xlim([Data.domain(1), Data.domain(2)])
                ylim([Data.domain(3), Data.domain(4)])
                
            else

                 %%  DG-solution scatter plot
                if Comp == 1
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}(:,1), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{h})_x$'}),'Interpreter','latex','Fontsize',14);
                else
                    scatter(Xh.Solution{1}, Xh.Solution{2}, 15, Xh.Solution{i}(:,2), 'filled');
                    title(strcat({'$('},{Xh.StrPlot{i}},{'^\mathrm{h})_y$'}),'Interpreter','latex','Fontsize',14);
                end
                colorbar
                hold on
                if Data.PlotGridSol
                    for kk = 1:region.ne
                        plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k')
                    end
                end
                axis equal
                xlim([Data.domain(1), Data.domain(2)])
                ylim([Data.domain(3), Data.domain(4)])

            end
        end
    end

    pause(1e-3)

end

end
