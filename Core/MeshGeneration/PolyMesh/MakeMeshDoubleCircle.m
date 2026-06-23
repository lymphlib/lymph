%> @file  MakeMeshDoubleCircle.m
%> @author Ilario Mazzieri
%> @date 5 June 2026
%> @brief Construction of a polygonal mesh for a circular domain.
%>  
%> Creation of a polygonal mesh for a circular domain
%> \cite quarteroni2009numerical (Example 5.1)
%>
%>   IDs for volume and boundary elements
%>           ___    ___
%>          /   \  /   \
%>         /     \/     \ 
%>     3  (      1       ) 2
%          \     /\     /
%>          \___/  \___/     
%>                
%>
%==========================================================================
%> @section classMakeMeshDoubleCircle Class description
%==========================================================================
%> @brief  Construction of the polygonal mesh (uses polymesher functions
%>contained in the Polymesher folder 
%>
%> @param Data                    Struct with problem's data
%> @param N                       Number of mesh elements
%> @param DomainLimits            Domain limits
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%>
%==========================================================================
function [FileNameOut] = MakeMeshDoubleCircle(Data,N,DomainLimits,FolderName,FileName) 

%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati

%xmin xmax ymin ymax
Dati.domain = DomainLimits;
Data.domain = DomainLimits;
Dati.Circle1.Radius = 0.75;
Dati.Circle1.Center = [-0.5 0];
Dati.Circle2.Radius = 0.75;
Dati.Circle2.Center = [0.5 0];


% polygonal mesh
[region] = PolyMesher(@DoubleCircle,N,500);    



%% Compute the neighbor structure
[region,neighbor] = MakeNeighbor(Data,region);

%% Otuput 
FileNameOut = [FolderName,'/',FileName];
save(FileNameOut,'region','neighbor');