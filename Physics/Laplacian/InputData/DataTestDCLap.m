%% Poisson problem with Dirichlet conditions

Data.name = 'ConvTestDCLap';

Data.TagElLap = 1;         % Element tag
Data.TagBcLap = [2 3]; % Boundary tag
Data.LabBcLap = 'DD';    % Dirichlet/Neumann/Abso

%% Geometrical properties 
Data.domain       = [-1.5 1.5 -1 1]; % domain bounds for a new mesh
Data.N            = 1000;        % number of elements for a new mesh
Data.MeshFromFile = true;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'DataTestDCLap_1000_el'; %filename for mesh 

%% Space discretization
Data.degree  = 1;   % Polynomial degree
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
Data.u_ex     =  {@(x,y) 0.*x.*y};

% Gradient of the Exact Solution
Data.gradu_ex =  {@(x,y)  0.*x.*y; ...
                  @(x,y)  0.*x.*y};

