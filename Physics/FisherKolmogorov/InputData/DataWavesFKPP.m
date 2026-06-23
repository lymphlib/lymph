%% Convergence test for FK equation

Data.name = 'WavesTest';

Data.LabEl = {'K'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'NNNN'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 0;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain          = [-1 3 0 1];
Data.N               = 1000; 
Data.meshfromfile    = false;
Data.foldername      = 'InputMesh';
Data.VTKMeshFileName = 'WavesMesh.vtk';
Data.meshfileseq     = "FK";

Data.isotropy    = true;
Data.AxnDiffFile = "...";

%% Material properties 

alpha           = 1;
d               = 1e-3;

Data.D_ext       = {@(x,y,t) d+0.*x};  % External diffusion
Data.D_axn       = {@(x,y,t) 0.*x};       % Axonal diffusion
Data.alpha       = {@(x,y,t) alpha+0.*x};     % Reaction parameter

% Forcing Term
Data.homog_source_f = true;
Data.source_f       = {@(x,y,t,D,alpha) 0.*x};

% Boundary Conditions
Data.DirBC    =  {@(x,y,t) 0.25*(1+tanh(8+5/12*alpha.*t-sqrt(alpha/(24*d)).*x)).^2};

% Exact Solution (if any)
Data.c_ex     =  {@(x,y,t) 0.25*(1+tanh(8+5/12*alpha.*t-sqrt(alpha/(24*d)).*x)).^2};

% Gradient of the Exact Solution
Data.dc_dx_ex =  {@(x,y,t) -0.5*sqrt(alpha/(24*d))*(sech(8+5/12*alpha.*t-sqrt(alpha/(24*d)).*x)).^2.*(1+tanh(8+5/12*alpha.*t-sqrt(alpha/(24*d)).*x))};
Data.dc_dy_ex =  {@(x,y,t) 0.*x};
Data.dc_dt_ex =  {@(x,y,t) 5/24*alpha*(sech(8+5/12*alpha.*t-sqrt(alpha/(24*d)).*x)).^2.*(1+tanh(8+5/12*alpha.*t-sqrt(alpha/(24*d)).*x))};

%% Discretization properties

%% Time discretization
Data.t0    = 0;
Data.T     = 1.00;
Data.dt    = 1e-2;

Data.theta    = 0.5;

%% Space discretization
Data.degree        = 3;		        % Polynomial degree
Data.penalty_coeff = 10;   			% Penalty coefficient

%% Nonlinear solver settings
Data.NonLinearSolver = 'Picard iterations';     % Nonlinearity treatment: Semi-implicit or Picard iterations

Data.NLS_tolerance  = 1e-10;
Data.NLS_max_it     = 1000;
Data.NLS_relax      = 1;

%% Visualization settings
Data.PlotExact          = true;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 40;
Data.NPtsVisualization  = 4;

%% Assembling settings
Data.quadrature = 'QF';           % Method of integrals calculation: Quadrature-Free - Subtriangulation

%% Adaptivity
Data.Adaptivity         = true;                % Flag of adaptivity                                             
Data.maxDegree          = 3;                    % Maximum polynomial degree
Data.AdaptFunc          = @(tau_r) floor(1 + (Data.maxDegree-1)*tanh(5*tau_r));     % Adaptivity function
Data.AdaptIts           = 2;                   % Maximum iterations of adaptivity
Data.AdaptivityStep     = 10; 
