%% Convergence test for Heat equation

Data.name = 'hConvergenceTest';

Data.TagElLap = 1;         % Element tag
Data.TagBcLap = [2 3 4 5]; % Boundary tag
Data.LabBcLap = 'DDDD';    % Dirichlet/Neumann/Abso

Data.TagApplyBCs = 1;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain             = [0 1 0 1];
Data.N                  = 1000;
Data.MeshFromFile       = true;
Data.FolderName         = 'InputMesh';
Data.VTKMeshFileName    = 'Mesh.vtk';
Data.meshfileseq        = 'UnitSquare_30_el.mat';

%% Material properties 

Data.mu          = {@(x,y,t) 1+0.*x};     % Diffusion parameter
Data.sigma       = {@(x,y,t) 0+0.*x};     % Reaction parameter

% Forcing Term
Data.homog_source_f = false;
Data.source_f       = {@(x,y,t) -(cos(pi*x).*cos(pi*y)+2).*exp(-t)+2*pi*pi*cos(pi*x).*cos(pi*y).*exp(-t)};

% Boundary Conditions
Data.DirBC    = @(x,y,t) (cos(pi*x).*cos(pi*y)+2).*exp(-t);

% Exact Solution (if any)
Data.u_ex     =  @(x,y,t) (cos(pi*x).*cos(pi*y)+2).*exp(-t);

% Gradient of the Exact Solution
Data.du_dx_ex =  @(x,y,t) -pi*sin(pi*x).*cos(pi*y).*exp(-t);
Data.du_dy_ex =  @(x,y,t) -pi*cos(pi*x).*sin(pi*y).*exp(-t);
Data.du_dt_ex =  @(x,y,t) -(cos(pi*x).*cos(pi*y)+2).*exp(-t);

%% Discretization properties

%% Time discretization

Data.t0     = 0;
Data.T      = 1e-4;
Data.dt     = 1e-5;
Data.theta  = 0.5;

%% Space discretization

Data.degree        = [1, 2, 3, 4, 5, 6, 7, 8];             % Polynomial degree
Data.penalty_coeff = 10;                                   % Penalty coefficient

%% Visualization settings

Data.PlotExact          = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 2e-3;

%% Save solution settings

Data.SaveSolutionStep = 0.5;
