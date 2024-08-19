%> @file  WriteVtk.m
%> @author Ilario Mazzieri, Mattia Corti, Stefano Bonetti
%> @date 24 July 2024
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
%> @param append     Boolean variable
%>                   if 'true', append to an existing .vtk file
%>                   if 'false', create a new .vtk file and write on it
%> 
%> @retval []
%>
%==========================================================================

function WriteVtk(fname, x, y, z, val, conn, prop_name, append)

% Check if we save a scalar- or vector- field
n_val = size(val,2);

if append

    fid = fopen(fname, 'a');

    if n_val == 1
        % Scalar value
        fprintf(fid,['SCALARS  ', prop_name,' float 1\n']);
        fprintf(fid,'LOOKUP_TABLE default\n');
        fprintf(fid,'%f\n',val);
        
    else        
        % Vectors
        fprintf(fid,['VECTORS  ' ,prop_name, '  float\n']);
        val = [val zeros(size(val,1),1)];
        fprintf(fid,'%f %f %f\n',val');
    end

else

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
    coords = [x y z];
    fprintf(fid,'%f %f %f\n',coords');
    
    numpoly = size(conn,1);
    numval = size(conn,1)*size(conn,2) + numpoly;
    
    % Print of the grid connectivity
    fprintf(fid,'CELLS %i %i\n',numpoly, numval);
    connout = [3*ones(size(conn,1),1) conn-1];
    fprintf(fid,'%i %i %i %i\n',connout');
    
    % Cell types: 7 is "Polygon"
    fprintf(fid,'CELL_TYPES %i\n', numpoly);
    out = repmat('7\n',1,numpoly);
    fprintf(fid,out);
    
    % Point data
    fprintf(fid,'POINT_DATA %i\n', length(x));
    
    if n_val == 1
        % Scalar value
        fprintf(fid,['SCALARS  ', prop_name,' float 1\n']);
        fprintf(fid,'LOOKUP_TABLE default\n');
        fprintf(fid,'%f\n',val);
        
    else        
        % Vectors
        fprintf(fid,['VECTORS  ' ,prop_name, '  float\n']);
        val = [val zeros(size(val,1),1)];
        fprintf(fid,'%f %f %f\n',val');
    end

end

fclose(fid);