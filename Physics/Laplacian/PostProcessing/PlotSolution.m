%> @file  PlotSolution.m
%> @author Ilario Mazzieri
%> @date 16 April 2023
%> @brief Plot the solutions (approximated and exact).
%>
%==========================================================================
%> @section classPlotSolution Class description
%==========================================================================
%> @brief Simple plot of the solutions
%>
%> @param Data        Struct with problem's data
%> @param UPlot  Solution array having the following structure
%>               Uplot = [ x | y | uh(x,y) | uex(x,y)];
%>
%> @retval ~ 
%>
%==========================================================================
function [] = PlotSolution(Data,UPlot)

if(~isempty(Data.TagElLap))   
    figure;
    tri = delaunay(UPlot(:,1),UPlot(:,2));
    subplot(1,2,1)
    hfig1 = trisurf(tri, UPlot(:,1), UPlot(:,2), UPlot(:,3));
    title('u_h(x,y)');
    axis vis3d;
    lfig1 = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
    subplot(1,2,2)
    hfig2 = trisurf(tri, UPlot(:,1), UPlot(:,2), UPlot(:,4));
    title('u_{ex}(x,y)')
    axis vis3d;
    lfig2 = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
end