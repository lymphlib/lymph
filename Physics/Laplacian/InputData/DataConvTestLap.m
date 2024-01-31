%% Poisson problem in [0,1]^2 with Dirichlet conditions

Data.name = 'Lap';

Data.TagElLap = 1;         % Element tag
Data.TagBcLap = [2 3 4 5]; % Boundary tag
Data.LabBcLap = 'DDDD';    % Dirichlet/Neumann/Abso

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = {100, 200, 400, 800};        % number of elements for a new mesh
Data.MeshFromFile = true;     % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = {"Lap_100_el.mat","Lap_200_el.mat", ... 
                     "Lap_400_el.mat", "Lap_800_el.mat"};  %filename for mesh 

%% Space discretization
Data.degree  = 4;   % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Visualization settings
Data.PlotExact   = true;
Data.PlotGridSol = true;

%% Material properties 
Data.mu       = {@(x,y) 1.*x.^0.*y.^0};

% Forcing Term
Data.source   = {@(x,y) 2*(4*pi^2)*sin(2*pi*x).*cos(2*pi*y)};

% Boundary Conditions
Data.DirBC    = {@(x,y) sin(2*pi*x).*cos(2*pi*y)};

% Exact Solution (if any)
Data.u_ex     =  {@(x,y) sin(2*pi*x).*cos(2*pi*y)};

% Gradient of the Exact Solution
Data.gradu_ex =  {@(x,y)  2*pi*cos(2*pi*x).*cos(2*pi*y); ...
                  @(x,y) -2*pi*sin(2*pi*x).*sin(2*pi*y)};

