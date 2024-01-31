%> @file  MakeMeshQuadridomain.m
%> @author Ilario Mazzieri
%> @date 8 March 2023
%> @brief Construction of a polygonal mesh for a rectangular domain.
%>
%> Creation of a polygonal mesh for a rectangular domain
%>  Omega = [xmin, xmax] x [ymin, ymax] divided into four rectangles
%>  sharing the interfaces {0}x[ymin,ymax] and [xmin xmax]x{0}
%>
%>           IDs for volume and boundary elements
%>            _10__ ___9__
%>           |     |      |
%>        11 |  2  |   3  |8
%>           |_____|______|
%>           |     |      |
%>        12 |  1  |   4  |7
%>           |_____|______|
%>              5      6
%>
%======================================================================
%> @section classMakeMeshQuadridomain Class description
%======================================================================
%> @brief  Construction of the polygonal mesh (uses polymesher functions
%>contained in the Polymesher folder
%>
%> @param Data			  Data structure
%> @param N                       Number of mesh elements
%> @param DomainLimits		  Domain limits
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%> @param MeshType                String 'C' for cartesian grid, 'P' for
%polygonal grid
%> @param SimType		  Simulation type
%>
%> @retval FileNameOut            File name of the *.mat structure
%>                                containing mesh info
%>
%======================================================================
function [FileNameOut] = MakeMeshQuadridomain(Data,N,DomainLimits,FolderName,FileName,MeshType,SimType)

%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% Rectangular domain
global Dati

%xmin<0 xmax>0 ymin<0 ymax>0
Dati.domain = DomainLimits;
Data.domain = DomainLimits;


if (strcmp(MeshType,'P')==1)
    
    % polygonal mesh
    [region] = PolyMesher_QuadriDomain(@RectangleQD,N,100);
    
elseif (strcmp(MeshType,'C')==1)
    
    % cartesian mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    N = nelx*nely;
    
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax+dx/2:dx:bx, ay+dy/2:dy:by);
    P = [X(:) Y(:)];
    
    [region] = PolyMesher_QuadriDomainCartesian(@RectangleQD,N,100,P);
    
end

%% Compute the neighbor structure
[region,neighbor] = MakeNeighborQuadridomain(Data,region,SimType);

%% Otuput 
FileNameOut = [FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat'];
save(FileNameOut,'region','neighbor');
