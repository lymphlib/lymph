%% Convergence test for FK equation - h convergence

Data.name = 'ConvTestDir';


Data.LabEl = {'K'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'DDDD'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 1;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain             = [0 1 0 1];
Data.N                  = {100, 300, 800, 1000};
Data.meshfromfile       = false;
Data.foldername         = 'InputMesh';
Data.VTKMeshFileName    = 'Mesh.vtk';
Data.meshfileseq        = ["FK_100_el.mat", "FK_300_el.mat", "FK_800_el.mat", "FK_1000_el.mat"];

Data.isotropy           = true;
Data.AxnDiffFile        = "...";

%% Material properties 

Data.D_ext       = {@(x,y,t) 1+0.*x};     % External diffusion
Data.D_axn       = {@(x,y,t) 0+0.*x};     % Axonal diffusion
Data.alpha       = {@(x,y,t) 1+0.*x};   % Reaction parameter

% Forcing Term
Data.homog_source_f = false;
Data.source_f       = {@(x,y,t,D,alpha) -(cos(pi*x).*cos(pi*y)+2).*exp(-t)+2*pi*pi*D.*cos(pi*x).*cos(pi*y).*exp(-t)-alpha.*(cos(pi*x).*cos(pi*y)+2).*exp(-t).*(1-(cos(pi*x).*cos(pi*y)+2).*exp(-t))};

% Boundary Conditions
Data.DirBC     = {@(x,y,t) (cos(pi*x).*cos(pi*y)+2).*exp(-t)};
Data.gradDirBC = {@(x,y,t) -pi*sin(pi*x).*cos(pi*y).*exp(-t), @(x,y,t) -pi*cos(pi*x).*sin(pi*y).*exp(-t)};

% Exact Solution (if any)
Data.c_ex     = {@(x,y,t) (cos(pi*x).*cos(pi*y)+2).*exp(-t)};

% Gradient of the Exact Solution
Data.dc_dx_ex =  {@(x,y,t) -pi*sin(pi*x).*cos(pi*y).*exp(-t)};
Data.dc_dy_ex =  {@(x,y,t) -pi*cos(pi*x).*sin(pi*y).*exp(-t)};
Data.dc_dt_ex =  {@(x,y,t) -(cos(pi*x).*cos(pi*y)+2).*exp(-t)};

%% Discretization properties

%% Time discretization
Data.t0 = 0;
Data.T  = 1e-4;
Data.dt = 1e-5;

Data.theta   = 0.5;

%% Space discretization
Data.degree        = 2;             % Polynomial degree
Data.penalty_coeff = 10;            % Penalty coefficient

%% Nonlinear solver settings
Data.NonLinearSolver = 'Semi-implicit';     % Nonlinearity treatment: Semi-implicit or Picard iterations

Data.NLS_tolerance  = 1e-10;
Data.NLS_max_it     = 10;
Data.NLS_relax      = 1;

%% Visualization settings
Data.PlotExact          = true;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 1;
Data.NPtsVisualization  = 4;

%% Assembling settings
Data.quadrature = 'QF';           % Method of integrals calculation: Quadrature-Free - Subtriangulation

%% Adaptivity
Data.Adaptivity         = false;                % Flag of adaptivity                                             
Data.maxDegree          = 2;                    % Maximum polynomial degree
Data.AdaptFunc          = @(tau_r) floor(1 + 2*Data.maxDegree/pi* atan(tau_r));     % Adaptivity function
Data.AdaptIts           = 1;                   % Maximum iterations of adaptivity
Data.AdaptivityStep     = 1; 
