%> @file  MakeMeshBidomainVert.m
%> @author Ilario Mazzieri
%> @date 8 March 2023
%> @brief Construction of a polygonal mesh for a rectangular domain.
%>
%> Creation of a polygonal mesh for a rectangular domain
%>  Omega = [xmin, xmax] x [ymin, ymax] divided into two rectangles
%>           sharing the interface {0}x[ymin,ymax]
%>
%>           IDs for volume and boundary elements
%>            __7__ ___6__
%>           |     |      |
%>          8|  1  |   2  |5
%>           |_____|______|
%>              3      4
%>
%======================================================================
%> @section classMakeMeshBidomainVert Class description
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
function [FileNameOut] = MakeMeshBidomainVert(Data,N,DomainLimits,FolderName,FileName,MeshType,SimType) 

%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% Rectangular domain
global Dati

% ATTENTION!!! First y and then x
% ymin ymax xmin<0 xmax>0
Dati.domain = [DomainLimits(3) DomainLimits(4) DomainLimits(1) DomainLimits(2)];
Data.domain = DomainLimits;

if (strcmp(MeshType,'P')==1)
    
    % polygonal mesh
    [region] = PolyMesher_BiDomain(@RectangleBD,N,100);
    
elseif (strcmp(MeshType,'C')==1)
    
    % cartesian mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    N = nelx*nely;
    
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax+dx/2:dx:bx, ay+dy/2:dy:by);
    P = [X(:) Y(:)];
    
    [region] = PolyMesher_BiDomainCartesian(@RectangleBD,N,100,P);
    
end


%% Compute the neighbor structure
[region,neighbor] = MakeNeighborBidomainVert(Data,region,SimType);

%% Otuput 
FileNameOut = [FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat'];
save(FileNameOut,'region','neighbor');
