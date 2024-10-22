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
%> @param Data                    Structure with problem data
%> @param N                       Number of mesh elements
%> @param DomainLimits            Domain limits
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%> @param MeshType                String 'C' for cartesian grid, 'P' for
%polygonal grid
%> @param SimType                 Simulation type
%>
%> @retval FileNameOut            File name of the *.mat structure
%>                                containing mesh info
%>
%======================================================================
function [FileNameOut] = MakeMeshBidomainHor(Data,N,DomainLimits,FolderName,FileName,MeshType,SimType)
%% SET DIRECTORIES AND NAMES

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati

if nargin < 7
    SimType = 'laplacian';
    Data.TagBcLap = [3, 4, 5, 6, 7, 8];
    Data.LabBcLap = 'DDDDDD';
end


Dati.domain = DomainLimits;
if (strcmp(MeshType,'P')==1)
    
    % polygonal mesh
    [region] = PolyMesher_BiDomainHor(@Rectangle,N,100);
    
elseif (strcmp(MeshType,'C')==1)
    
    % cartesian mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    N = nelx*nely;
    
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax+dx/2:dx:bx, ay+dy/2:dy:by);
    P = [X(:) Y(:)];
    
    [region] = PolyMesher_BiDomainHorCartesian(@Rectangle,N,100,P);
    
end

%% Compute the neighbor structure
[region,neighbor] = MakeNeighborBidomainHor(Data,region,SimType);

%% Otuput 
FileNameOut = [FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat'];
save(FileNameOut,'region','neighbor');

