%% Poisson problem in [0,1]^2 with Dirichlet conditions

Data.name = 'TestLapDir';

Data.LabEl = {'L'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'DDDD'};  % (D)irichlet/(N)eumann/(A)bsorbing

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 100;      % number of elements for a new mesh
Data.MeshFromFile = false;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'Lap'; %filename for mesh

%% Space discretization
Data.degree  = 3;        % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Quadrature settings
Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = true;
Data.PlotGridSol       = true;
Data.NPtsVisualization = 20;

%% Material properties 
Data.mu       = {@(x,y) 1.*x.^0.*y.^0};

% Forcing Term
Data.source   = {@(x,y) 10*(4*pi^2)*sin(2*pi*x).*cos(2*pi*y)};

% Boundary Conditions
Data.DirBC    = {@(x,y) 5*sin(2*pi*x).*cos(2*pi*y)};

% Gradient of the Exact Solution
Data.gradDirBC =  {@(x,y)  10*pi*cos(2*pi*x).*cos(2*pi*y); ...
                  @(x,y) -10*pi*sin(2*pi*x).*sin(2*pi*y)};

% Exact Solution (if any)
Data.u_ex     =  {@(x,y) 5*sin(2*pi*x).*cos(2*pi*y)};

% Gradient of the Exact Solution
Data.gradu_ex =  {@(x,y)  10*pi*cos(2*pi*x).*cos(2*pi*y); ...
                  @(x,y) -10*pi*sin(2*pi*x).*sin(2*pi*y)};


%% Adaptivity
Data.Adaptivity = false;                 % Flag of adaptivity                                             
Data.maxDegree  = 3;                    % Maximum polynomial degree
Data.AdaptFunc  = @(tau_r) floor(1 + 2*Data.maxDegree/pi* atan(tau_r));     % Adaptivity function
Data.AdaptIts   = 5;                    % Maximum iterations of adaptivity
