%> @file  ScatterPlotIndicator.m
%> @author Mattia Corti, Caterina B Leimer Saglio
%> @date 10 February 2026
%> @brief Scatter visualization of the p-adaptivity indicator
%>
%==========================================================================
%> @section classScatterPlotIndicator Class description
%==========================================================================
%> @brief            Scatter visualization of the p-adaptivity indicator
%
%> @param Xh        Cell array that contains the approximate solutions. It has the following form:
%> Xh = { x , y , u_h (x,y), v_h(x,y), ... } for scalar-valued outputs
%> Xh = { x , y , [(u_h)_x(x,y) (u_h)_y(x,y)], [(v_h)_x(x,y) (v_h)_y(x,y)],
%... } for vector-valued outputs
%> @param region     Mesh region struct
%> @param Data       Struct with problem's data
%
%> @retval ~
%>
%==========================================================================

function ScatterPlotIndicator(Xh, region, Data)

        figure

        subplot(1,2,1);
        %%  scatter polynomial degree
        scatter(Xh.Solution{1}, Xh.Solution{2}, 1, Xh.Solution{strcmp(Xh.StrPlot,'p')}, 'filled');
        title(strcat({'$'},{'p_K'},{'^\mathrm{h}$'}),'Interpreter','latex');
        colorbar
        colormap("jet")
        hold on
        if Data.PlotGridSol
            for kk = 1:region.ne
                plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k', "LineWidth", 1)
            end
        end
        axis equal
        xlim([Data.domain(1), Data.domain(2)])
        ylim([Data.domain(3), Data.domain(4)])


        subplot(1,2,2);
        %%  scatter complete tau indicator
        scatter(Xh.Solution{1}, Xh.Solution{2}, 1, Xh.Solution{strcmp(Xh.StrPlot,'tau')}, 'filled');
        title(strcat({'$'},{'\tau_K'},{'^\mathrm{h}$'}),'Interpreter','latex');
        colorbar  
        colormap("jet")
        hold on
        if Data.PlotGridSol
            for kk = 1:region.ne
                plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k', "LineWidth", 1)
            end
        end
        axis equal
        xlim([Data.domain(1), Data.domain(2)])
        ylim([Data.domain(3), Data.domain(4)])

        pause(1e-3)

        %%  scatter different parts of tau
        figure

        subplot(1,4,1);
        scatter(Xh.Solution{1}, Xh.Solution{2}, 1, (Xh.Solution{strcmp(Xh.StrPlot,'tauE')}), 'filled');
        title(strcat({'$'},{'\tau_{K,e}'},{'^\mathrm{h}$'}),'Interpreter','latex');
        colorbar
        colormap("jet")
        hold on
        if Data.PlotGridSol
            for kk = 1:region.ne
                plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k', "LineWidth", 1)
            end
        end
        axis equal
        xlim([Data.domain(1), Data.domain(2)])
        ylim([Data.domain(3), Data.domain(4)])
        

        subplot(1,4,2);
        scatter(Xh.Solution{1}, Xh.Solution{2}, 1, Xh.Solution{strcmp(Xh.StrPlot,'tauJ')}, 'filled');
        title(strcat({'$'},{'\tau_{K,j}'},{'^\mathrm{h}$'}),'Interpreter','latex');
        colorbar
        colormap("jet")
        hold on
        if Data.PlotGridSol
            for kk = 1:region.ne
                plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k', "LineWidth", 1)
            end
        end
        axis equal
        xlim([Data.domain(1), Data.domain(2)])
        ylim([Data.domain(3), Data.domain(4)])
        

        subplot(1,4,3);
        scatter(Xh.Solution{1}, Xh.Solution{2}, 1, Xh.Solution{strcmp(Xh.StrPlot,'tauN')}, 'filled');
        title(strcat({'$'},{'\tau_{K,n}'},{'^\mathrm{h}$'}),'Interpreter','latex');
        colorbar
        colormap("jet")
        hold on
        if Data.PlotGridSol
            for kk = 1:region.ne
                plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k', "LineWidth", 1)
            end
        end
        axis equal
        xlim([Data.domain(1), Data.domain(2)])
        ylim([Data.domain(3), Data.domain(4)])
        

        subplot(1,4,4);
        scatter(Xh.Solution{1}, Xh.Solution{2}, 1, Xh.Solution{strcmp(Xh.StrPlot,'tauT')}, 'filled');
        title(strcat({'$'},{'\tau_{K,t}'},{'^\mathrm{h}$'}),'Interpreter','latex');
        colorbar
        colormap("jet")
        hold on
        if Data.PlotGridSol
            for kk = 1:region.ne
                plot([region.coords_element{kk}(:,1); region.coords_element{kk}(1,1)], [region.coords_element{kk}(:,2); region.coords_element{kk}(1,2)], 'k', "LineWidth", 1)
            end
        end
        axis equal
        xlim([Data.domain(1), Data.domain(2)])
        ylim([Data.domain(3), Data.domain(4)])
        
        pause(1e-3)
    end

   
   
