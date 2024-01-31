%> @file  MakeMeshMonodomain.m
%> @author Ilario Mazzieri
%> @date 8 March 2023 
%> @brief Construction of a polygonal mesh for a rectangular domain.
%>  
%> Creation of a polygonal mesh for a rectangular domain
%>   Omega = [xmin, xmax] x [ymin, ymax]
%>
%>   IDs for volume and boundary elements
%>            _____4_____
%>           |           |
%>          5|     1     |3   
%>           |___________|    
%>                 2
%>
%==========================================================================
%> @section classMakeMeshMonodomain Class description
%==========================================================================
%> @brief  Construction of the polygonal mesh (uses polymesher functions
%>contained in the Polymesher folder 
%>
%> @param Data                    Struct with problem's data
%> @param N                       Number of mesh elements
%> @param DomainLimits            Domain limits
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%> @param MeshType                String 'C' for cartesian grid, 'P' for
%polygonal grid
%> @param SimType                 String simulation type, used for boundary tag       
%>
%> @retval FileNameOut            File name of the *.mat structure 
%>                                containing mesh info 
%>
%==========================================================================
function [FileNameOut] = MakeMeshMonodomain(Data,N,DomainLimits,FolderName,FileName,MeshType,SimType) 

%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati

%xmin xmax ymin ymax
Dati.domain = DomainLimits;
Data.domain = DomainLimits;


if (strcmp(MeshType,'P')==1)
    
    % polygonal mesh
    [region] = PolyMesher(@Rectangle,N,100);    
    
elseif (strcmp(MeshType,'C')==1)
    
    % cartesian mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    N = nelx*nely;
    
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax+dx/2:dx:bx, ay+dy/2:dy:by);
    P = [X(:) Y(:)];
    
    [region] = PolyMesher(@Rectangle,N,100,P);
    
end



%% Compute the neighbor structure
[region,neighbor] = MakeNeighbor(Data,region,SimType);

%% Otuput 
FileNameOut = strcat(FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat');
save(FileNameOut,'region','neighbor');
