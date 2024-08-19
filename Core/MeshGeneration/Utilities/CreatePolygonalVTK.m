%> @file  CreatePolygonalVTK.m
%> @author Mattia Corti and Ilario Mazzieri
%> @date 4 May 2023
%> @brief Save the polygonal mesh as a polygonal VTK file 
%>
%==========================================================================
%> @section classCreatePolygonalVTK Class description
%==========================================================================
%> @brief            Save the polygonal mesh as a polygonal VTK file
%>
%> @param Data       Struct with problem's data
%> @param Setup      Struct with simulation's setup
%> @param region     Mesh region struct
%>
%==========================================================================

function CreatePolygonalVTK(Data, Setup, region)

    %% Creation of the VTK file
    VTKName = [Setup.OutFolderVTK, '/', Data.VTKMeshFileName];
    fileID = fopen(VTKName,'w');

    %% Date statement reconstruction
    DateState = ['surface written ',num2str(year(datetime)),'-',...
        num2str(month(datetime)),'-',day(datetime),'T', ...
        num2str(hour(datetime)),':',num2str(minute(datetime)),':',...
        num2str(floor(second(datetime))),'\n'];

    %% Creation of the VTK header file
    fprintf(fileID,"# vtk DataFile Version 2.0\n");
    fprintf(fileID,DateState);
    fprintf(fileID,"ASCII\n\n");
    fprintf(fileID,"DATASET POLYDATA\n");

    %% Points
    fprintf(fileID,"POINTS %d float\n",size(region.coord,1));
    app = [region.coord zeros(size(region.coord,1),1)];
    fprintf(fileID,'%d %6.7f %6.7f\r\n',app');

    s = region.ne+sum(cell2mat(cellfun(@length,region.connectivity,'UniformOutput',false)));

    %% Polygons
    fprintf(fileID,"\nPOLYGONS %d %d", region.ne, s);

    for ii = 1:length(region.connectivity)
        fprintf(fileID, "\n%d ", length(region.connectivity{ii}));
        fprintf(fileID, "%d ", region.connectivity{ii}-1);
    end

    fprintf(fileID,"\nCELL_DATA %d",region.ne);
    fprintf(fileID,"\nSCALARS POLYHEDRA int 1");
    fprintf(fileID,"\nLOOKUP_TABLE default ");
    fprintf(fileID,"%d ", region.id);
    

    %% Close the file
    fclose(fileID);
end