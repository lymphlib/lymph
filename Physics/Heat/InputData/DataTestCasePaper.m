%% Patch test for Heat equation

Data.name = 'TestCase';

Data.TagElLap = 1;         % Element tag
Data.TagBcLap = [2 3 4 5]; % Boundary tag
Data.LabBcLap = 'DDDD';    % Dirichlet/Neumann/Abso

Data.TagApplyBCs = 1;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain             = [-1.5 1.5 -1 1];
Data.N                  = 250;                     % Number of mesh elements
Data.MeshFromFile       = true;
Data.FolderName         = 'InputMesh';
Data.VTKMeshFileName    = 'Mesh.vtk';
Data.meshfileseq        = 'DataTestDCLap_250_el.mat'; % Names of mesh files

%% Material properties 

Data.mu          = {@(x,y) 0.1+0.*x};     % Diffusion parameter
Data.sigma       = {@(x,y) 0+0.*x};       % Reaction parameter

% Forcing Term
Data.homog_source_f = true;
Data.source_f       = {@(x,y,t) 0.*x};

% Boundary Conditions
Data.DirBC    =  @(x,y,t) (x>=0);

% Exact Solution (if any)
Data.u_ex     =  @(x,y,t) (x>=0);

% Gradient of the Exact Solution
Data.du_dx_ex =  @(x,y,t) 0.*x;
Data.du_dy_ex =  @(x,y,t) 0.*x;
Data.du_dt_ex =  @(x,y,t) 0.*x;

%% Discretization properties

%% Time discretization

Data.t0     = 0;
Data.T      = 1;
Data.dt     = 2e-3;
Data.theta  = 0.5;


%% Space discretization

Data.degree        = 5;          % Polynomial degree
Data.penalty_coeff = 10;         % Penalty coefficient

%% Quadrature settings

Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings

Data.PlotExact          = true;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 50;
Data.NqnVisualization   = 5;        % Data.NqnVisualization must be odd and strictly greater than 1 (see Gauleg.m) 


%% Save solution settings

Data.SaveSolutionStep = 1;
