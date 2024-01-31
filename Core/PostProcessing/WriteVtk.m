%> @file  WriteVtk.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 3 October 2023
%> @brief  Write VTK files
%>
%==========================================================================
%> @section classWriteVtk Class description
%> @brief  Write VTK files
%
%> @param fname      File name
%> @param x,y,z      Grid Coordinates
%> @param val        Value to be printed in output:
%>                   if 1 column, the name is prop_name,
%>                   if 2 columns, the names are prop_name_x, prop_name_y
%> @param conn       Connectivity matrix
%> @param prop_name  Property name
%>
%> @retval []
%>
%==========================================================================

function WriteVtk(fname, x, y, z, val, conn, prop_name)

n_val = size(val,2);

fid = fopen(fname, 'w');

% File header (requested by the file format)
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment\n');

% Type of file: ASCII
fprintf(fid, 'ASCII\n');

% Data type: in this case we created an unstructured grid
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Mesh points
fprintf(fid,'POINTS %i float\n', length(x));

for i=1:length(x)
    fprintf(fid,'%g %g %g \n', x(i), y(i), z(i));
end

numpoly = length(conn);
numval = 0;

for i=1:numpoly
    numval = numval + length(conn(i,:)) + 1;
end

% Print of the grid connectivity
fprintf(fid,'CELLS %i %i\n',numpoly, numval);
for i=1:numpoly
    numpoint = length(conn(i,:));
    fprintf(fid,'%i ',[numpoint , conn(i,:)-1']);
    fprintf(fid,'\n');
end

% Cell types: 7 is "Polygon"
fprintf(fid,'CELL_TYPES %i\n', numpoly);
for i=1:numpoly
    fprintf(fid,'7\n' );
end

% Point data
fprintf(fid,'POINT_DATA %i\n', length(x));

if n_val == 1
    % Scalar value
    fprintf(fid,['SCALARS  ', prop_name,' float 1\n']);
    fprintf(fid,'LOOKUP_TABLE default\n');

    for i=1:length(x)
        fprintf(fid,'%g\n', val(i));
    end
    
else

    %scalars
    fprintf(fid,['SCALARS  ', prop_name,'_x float 1\n']);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:length(x)
        fprintf(fid,'%g\n', val(i,1));
    end
    fprintf(fid,['SCALARS  ', prop_name,'_y float 1\n']);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:length(x)
        fprintf(fid,'%g\n', val(i,2));
    end

    %vectors
    fprintf(fid,['VECTORS  ' ,prop_name, '  float\n']);
    for i=1:length(x)
        fprintf(fid,'%g %g %g\n',val(i,1), val(i,2), 0);
    end
end

fclose(fid);
