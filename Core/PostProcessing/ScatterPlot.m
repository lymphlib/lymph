%> @file  ScatterPlot.m
%> @author Mattia Corti and Ilario MAzzieri
%> @date 3 May 2023
%> @brief Scatter visualization of the elastic solution
%>
%==========================================================================
%> @section classScatterPlot Class description
%==========================================================================
%> @brief            Scatter visualization of the elastic solution
%
%> @param Gu         Matrix of either one of the following columns:
%> Gu = [ x | y |  u_h (x,y)   |  u_ex (x,y) ]
%> Gu = [ x | y | (u_h)_x(x,y) | (u_h)_y(x,y) | (u_ex)_x(x,y) | (u_ex)_y(x,y) ]
%> where x, y are the coordinates of the quadrature nodes,
%>          u_h is the numerical solution (either scalar or vectorial) in those nodes,
%>          u_ex is the exact solution (either scalar or vectorial)
%> @param Data       Struct with problem's data
%> @param region     Mesh region struct
%> @param Tag        Name to identify displ or velocity
%> @param Comp       Componenent of the solution (x or y) to be printed
%> (only in the vectorial case)
%
%> @retval ~
%>
%==========================================================================

function ScatterPlot(Gu, region, Data, Tag, Comp)

%% Check if data is scalar or vectorial
if ~(size(Gu,2)==4 || size(Gu,2)==6)
    error(strcat('ScatterPlot: first argument must have either 4 or 6 columns (found ', size(Gu,2),')'))
end
isScalar = (size(Gu,2) == 4);

figure

%% Control if plot the exact solution
if Data.PlotExact
    %% Exact-solution scatter plot
    subplot(1,3,2)
    if isScalar
        scatter(Gu(:,1), Gu(:,2), 15, Gu(:,4), 'filled');
        title(strcat({'$'},{Tag},{'^\mathrm{ex}$'}),'Interpreter','latex','Fontsize',14);
    else
        if Comp == 1
            scatter(Gu(:,1), Gu(:,2), 15, Gu(:,5), 'filled');
            title(strcat({'$('},{Tag},{'^\mathrm{ex})_x$'}),'Interpreter','latex','Fontsize',14);
        else
            scatter(Gu(:,1), Gu(:,2), 15, Gu(:,6), 'filled');
            title(strcat({'$('},{Tag},{'^\mathrm{ex})_y$'}),'Interpreter','latex','Fontsize',14);
        end
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
    if isScalar
        scatter(Gu(:,1), Gu(:,2), 15, Gu(:,3)-Gu(:,4), 'filled');
        title(strcat({'$'},{Tag},{'^\mathrm{h}-'},{Tag},{'^\mathrm{ex}$'}),'Interpreter','latex','Fontsize',14);
    else
        if Comp == 1
            scatter(Gu(:,1), Gu(:,2), 15, Gu(:,3)-Gu(:,5), 'filled');
            title(strcat({'$('},{Tag},{'^\mathrm{h})_x-('},{Tag},{'^\mathrm{ex})_x$'}),'Interpreter','latex','Fontsize',14);
        else
            scatter(Gu(:,1), Gu(:,2), 15, Gu(:,4)-Gu(:,6), 'filled');
            title(strcat({'$('},{Tag},{'^\mathrm{h})_y-('},{Tag},{'^\mathrm{ex})_y$'}),'Interpreter','latex','Fontsize',14);
        end
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

%%  DG-solution scatter plot
subplot(1,3,1)
if isScalar
    scatter(Gu(:,1), Gu(:,2), 15, Gu(:,3), 'filled');
    title(strcat({'$'},{Tag},{'^\mathrm{h}$'}),'Interpreter','latex','Fontsize',14);
else
    if Comp == 1
        scatter(Gu(:,1), Gu(:,2), 15, Gu(:,3), 'filled');
        title(strcat({'$('},{Tag},{'^\mathrm{h})_x$'}),'Interpreter','latex','Fontsize',14);
    else
        scatter(Gu(:,1), Gu(:,2), 15, Gu(:,4), 'filled');
        title(strcat({'$('},{Tag},{'^\mathrm{h})_y$'}),'Interpreter','latex','Fontsize',14);
    end
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

pause(1e-3)

end
