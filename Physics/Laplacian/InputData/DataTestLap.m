%% Poisson problem in [0,1]^2 with Dirichlet conditions

Data.name = 'ConvTestLapDir';

Data.TagElLap = 1;         % Element tag
Data.TagBcLap = [2 3 4 5]; % Boundary tag
Data.LabBcLap = 'DDDD';    % Dirichlet/Neumann/Abso

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 19600;     % number of elements for a new mesh
Data.MeshFromFile = true;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'Quad_19600_el.mat'; %filename for mesh

%% Space discretization
Data.degree  = 5;        % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Quadrature settings
Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = true;
Data.PlotGridSol       = true;
Data.NPtsVisualization = 5;

%% Material properties 
Data.mu       = {@(x,y) 1.*x.^0.*y.^0};

% Forcing Term
Data.source   = {@(x,y) 2*(400*pi^2)*sin(20*pi*x).*cos(20*pi*y)};

% Boundary Conditions
Data.DirBC    = {@(x,y) sin(20*pi*x).*cos(20*pi*y)};

% Exact Solution (if any)
Data.u_ex     =  {@(x,y) sin(20*pi*x).*cos(20*pi*y)};

% Gradient of the Exact Solution
Data.gradu_ex =  {@(x,y)  20*pi*cos(20*pi*x).*cos(20*pi*y); ...
                  @(x,y) -20*pi*sin(20*pi*x).*sin(20*pi*y)};

