%> @file  MakeMeshFlowAroundCilinder.m
%> @author Ilario Mazzieri
%> @date 8 October 2023 
%> @brief Construction of a polygonal mesh for a rectangular 
%>        domain with a circular inclusion.
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
%> @section classMakeMeshFlowAroundCilinder Class description
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
function [FileNameOut] = MakeMeshAcousticCavity(Data,N,DomainLimits,Circle,FolderName,FileName) 

%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati 

%xmin xmax ymin ymax
Dati.domain = DomainLimits;
Dati.circle = Circle;
Data.domain = DomainLimits;

% Rectangular domain in background
% P1 = PolyMshrRndPtSet(N,@AcousticBackground);
[~,P1] = PolyMesher(@AcousticBackground,N,100);    


% Structured points for circular domain
Dati.tol  = Dati.circle.Radius/4;
TetaPoints = 72;
rho = linspace(Dati.circle.Radius - Dati.tol, Dati.circle.Radius + Dati.tol,6);
p = [];
for i = 1 : length(rho)
    for teta = linspace(0,2*pi-2*pi/TetaPoints,TetaPoints)
        x = Dati.circle.Center(1) + rho(i)*cos(teta) ;
        y = Dati.circle.Center(2) + rho(i)*sin(teta) ;
        p = [p; x y];
    end
end

dcenter = sqrt((P1(:,1)-Dati.circle.Center(1)).^2 + (P1(:,2)-Dati.circle.Center(2)).^2);
P1 = P1(find(dcenter>(Dati.circle.Radius + Dati.tol)),:);

% Acoustic Cavity
%P2 = PolyMshrRndPtSet(50,@AcousticCavity);
[~,P2] = PolyMesher(@AcousticCavity,300,100);    


P = [P1;p;P2];

figure(1)
plot(P1(:,1),P1(:,2),'bo'); hold on;
plot(P2(:,1),P2(:,2),'ro'); hold on;
plot(p(:,1),p(:,2),'go'); hold on;

   
% polygonal mesh
figure;
[region] = PolyMesherAcousticCavity(@AcousticBackgroundAndCavity,size(P,1),1,P,Dati);    
% Plot the circular boundary

hold on;
p = [];
for teta = linspace(0,2*pi,TetaPoints)
        x = Dati.circle.Center(1) + Dati.circle.Radius*cos(teta) ;
        y = Dati.circle.Center(2) + Dati.circle.Radius*sin(teta) ;
        p = [p; x y];
end
plot(p(:,1),p(:,2),'g','Linewidth',2);

%% Compute the neighbor structure
[region,neighbor] = MakeNeighbor(Data,region,'waves');

%% Otuput 
FileNameOut = [FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat'];
save(FileNameOut,'region','neighbor');
