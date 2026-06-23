%% Patch test for Heat equation

Data.name = 'TestCase';

Data.LabEl = {'L'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'DDDD'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 1;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain             = [0 1 0 1];
Data.N                  = 250;                     % Number of mesh elements
Data.MeshFromFile       = false;
Data.FolderName         = 'InputMesh';
Data.VTKMeshFileName    = 'Mesh.vtk';
Data.meshfileseq        = 'DataTest'; % Names of mesh files

%% Material properties 

Data.mu          = {@(x,y,t) 1+0.*x};     % Diffusion parameter
Data.sigma       = {@(x,y,t) 0+0.*x};     % Reaction parameter

% Forcing Term
Data.homog_source_f = false;
Data.source_f       = {@(x,y,t) -(cos(pi*x).*cos(pi*y)+2).*exp(-t)+2*pi*pi*cos(pi*x).*cos(pi*y).*exp(-t)};

% Boundary Conditions
Data.DirBC        = {@(x,y,t) (cos(pi*x).*cos(pi*y)+2).*exp(-t)};
Data.gradDirBC    = {@(x,y,t) -pi.*(sin(pi*x).*cos(pi*y)).*exp(-t) ; @(x,y,t) -pi.*(cos(pi*x).*sin(pi*y)).*exp(-t)};

% Exact Solution (if any)
Data.u_ex     =  {@(x,y,t) (cos(pi*x).*cos(pi*y)+2).*exp(-t)};

% Gradient of the Exact Solution
Data.du_dx_ex =  {@(x,y,t) -pi*sin(pi*x).*cos(pi*y).*exp(-t)};
Data.du_dy_ex =  {@(x,y,t) -pi*cos(pi*x).*sin(pi*y).*exp(-t)};
Data.du_dt_ex =  {@(x,y,t) -(cos(pi*x).*cos(pi*y)+2).*exp(-t)};

%% Discretization properties

%% Time discretization

Data.t0     = 0;
Data.T      = 1;
Data.dt     = 2.5e-2;
Data.theta  = 0.5;

%% Space discretization

Data.degree        = 3;          % Polynomial degree
Data.penalty_coeff = 10;         % Penalty coefficient

%% Quadrature settings

Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings

Data.PlotExact          = true;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 10;
Data.NPtsVisualization  = 20; 

%% Adaptivity
Data.Adaptivity         = false;                 % Flag of adaptivity                                             
Data.maxDegree          = 5;                    % Maximum polynomial degree
Data.AdaptFunc          = @(tau_r) floor(1 + 2*Data.maxDegree/pi* atan(tau_r));     % Adaptivity function
Data.AdaptIts           = 2;                   % Maximum iterations of adaptivity
Data.AdaptivityStep     = 10; 