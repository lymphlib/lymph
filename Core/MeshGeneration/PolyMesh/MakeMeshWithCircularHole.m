%> @file  MakeMeshWithCircularHole.m
%> @author Ilario Mazzieri
%> @date 5 June 2026
%> @brief Construction of a polygonal mesh for a rectangular 
%>        domain with a circular hole.
%>  
%> Creation of a polygonal mesh for a rectangular domain
%>   Omega = [xmin, xmax] x [ymin, ymax] \ (x-xc)^2 + (y-yc)^2 <= r^2 
%>
%>   IDs for volume and boundary elements
%>            _____4_________
%>           |               |
%>          5|  o(6)   1     |3   
%>           |_______________|    
%>                 2
%>
%==========================================================================
%> @section classMakeMeshWithCircularHole Class description
%==========================================================================
%> @brief  Construction of the polygonal mesh (uses polymesher functions
%>contained in the Polymesher folder 
%>
%> @param Data                    Data structure
%> @param N                       Number of mesh elements
%> @param DomainLimits            Domain limits
%> @param Circle                  Struct with fields Radius and Center
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%>
%> @retval FileNameOut            File name of the *.mat structure 
%>                                containing mesh info 
%>
%==========================================================================
function [FileNameOut] = MakeMeshWithCircularHole(Data,N,DomainLimits,Circle,FolderName,FileName) 

%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati
Dati = struct;

%xmin xmax ymin ymax
Dati.domain = DomainLimits;
Dati.circle = Circle;
Data.domain = DomainLimits;

   
% polygonal mesh
[region] = PolyMesher(@CircleInclusion,N,100);    
    

%% Compute the neighbor structure
[region,neighbor] = MakeNeighbor(Data,region);

%% Otuput 
FileNameOut = [FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat'];
save(FileNameOut,'region','neighbor');
