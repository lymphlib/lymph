%> @file  MakeMeshBidomainHor.m
%> @author Ilario Mazzieri
%> @date 8 March 2023
%> @brief Construction of a polygonal mesh for a rectangular domain.
%>
%> Creation of a polygonal mesh for a rectangular domain
%>  \f$Omega = [xmin, xmax] \times [ymin, ymax]\f$ divided into two rectangles
%>  sharing the interface [xmin,xmax] x {0}.
%>
%>           IDs for volume and boundary elements
%>            __4__
%>           |     |
%>          5|  1  |3
%>           |_____|
%>           |     |
%>          6|  2  |8
%>           |_____|
%>              7
%>
%======================================================================
%> @section classMakeMeshBidomainHor Class description
%======================================================================
%> @brief  Construction of the polygonal mesh (uses polymesher functions
%>contained in the Polymesher folder
%>
%> @param N                       Number of mesh elements
%> @param DomainLimits            Domain limits
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%> @param MeshType                String 'C' for cartesian grid, 'P' for
%polygonal grid
%>
%> @retval FileNameOut            File name of the *.mat structure
%>                                containing mesh info
%>
%======================================================================
function [FileNameOut] = MakeMeshBidomainHor(N,DomainLimits,FolderName,FileName,MeshType)
%% SET DIRECTORIES AND NAMES

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati

Dati.domain = DomainLimits;
if (strcmp(MeshType,'P')==1)
    
    % polygonal mesh
    [region] = PolyMesher_BiDomainHor(@RectangleBD,N,100);
    
elseif (strcmp(MeshType,'C')==1)
    
    % cartesian mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    N = nelx*nely;
    
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax+dx/2:dx:bx, ay+dy/2:dy:by);
    P = [X(:) Y(:)];
    
    [region] = PolyMesher_BiDomainHorCartesian(@RectangleBD,N,100,P);
    
end

%% plot polygonal mesh and element numbers

PlotPolyMesh(region);

% hold on;
% for i = 1 : size(region.coord,1)
%     text(region.coord(i,1),region.coord(i,2), num2str(i));
% end

%% CONNECTIVITY - NEIGHBOURS AND BOUNDARY EDGES

[neighbour] = neighbours(region);

fid = fopen([FolderName,'B-connectivity'],'w');
k = 1;
for i = 1 : region.ne
    E_tag(i,1) = region.id(i);
    for j = 1 : neighbour.nedges(i)
        id_edge = neighbour.neigh{i}(j);
        if (E_tag(i,1) == 1 && id_edge == -1)
            %right edge
            B_tag = 3;
            if(j < neighbour.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
            %upper edge
            if(abs(region.coord(edge(1),2) - Dati.domain(4)) < 1.e-8 && ...
                    abs(region.coord(edge(2),2) - Dati.domain(4)) < 1.e-8)
                B_tag = 4;
            end
            %left edge
            if(abs(region.coord(edge(1),1) - Dati.domain(1)) < 1.e-8 && ...
                    abs(region.coord(edge(2),1) - Dati.domain(1)) < 1.e-8)
                B_tag = 5;
            end
            
            fprintf(fid,'%i  %i  %i\n', [B_tag, edge(1), edge(2)]);
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
        elseif (E_tag(i,1) == 2 && id_edge == -1)
            %left edge
            B_tag =6;
            if(j < neighbour.nedges(i))
                edge = [region.connectivity{i}(j) region.connectivity{i}(j+1)];
            else
                edge = [region.connectivity{i}(j) region.connectivity{i}(1)];
            end
            %right edge
            if(abs(region.coord(edge(1),1) - Dati.domain(2)) < 1.e-8 && ...
                    abs(region.coord(edge(2),1) - Dati.domain(2)) < 1.e-8)
                B_tag = 8;
            end
            
            %bottom edge
            if(abs(region.coord(edge(1),2) - Dati.domain(3)) < 1.e-8 && ...
                    abs(region.coord(edge(2),2) - Dati.domain(3)) < 1.e-8)
                B_tag = 7;
            end
            
            
            fprintf(fid,'%i  %i  %i\n', [B_tag, edge(1), edge(2)]);
            region.connectivity_bc(k,1:2) = [edge(1), edge(2)];
            region.bc_tag(k) = B_tag;
            k = k +1;
            
            
        end
        
        
    end
    
end
fclose(fid);

N_poly = length(E_tag);
Mesh.N     = N_poly;
Mesh.E_tag = E_tag;
Mesh.B_tag = B_tag;
Mesh.region = region;
Mesh.neighbor  = neighbour;
Mesh.domain_lim = Dati.domain;

FileNameOut = [FolderName,'/',FileName,'_',num2str(N_poly),'_el.mat'];

save(FileNameOut,'-struct','Mesh');
