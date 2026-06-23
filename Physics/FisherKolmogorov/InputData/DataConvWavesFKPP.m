%% Convergence test for FK equation - p convergence

Data.name = 'FisherKPP_p_convergence';

Data.LabEl = {'K'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'NNNN'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 0;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain          = [-1 3 0 1];
Data.N               = 300; 
Data.meshfromfile    = false;
Data.foldername      = 'InputMesh';
Data.VTKMeshFileName = 'WavesMesh.vtk';
Data.meshfileseq     = "FK";

Data.isotropy    = true;
Data.AxnDiffFile = "...";

%% Material properties 

alpha = 1;
d     = 1e-3;

Data.D_ext       = {@(x,y,t) d+0.*x};  % External diffusion
Data.D_axn       = {@(x,y,t) 0.*x};       % Axonal diffusion
Data.alpha       = {@(x,y,t) alpha+0.*x};     % Reaction parameter

% Forcing Term
Data.homog_source_f = true;
Data.source_f   = @(x,y,t,D,alpha) 0.*x;

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
Data.T     = 0.10;
Data.dt    = 1e-3;

Data.theta    = 0.5;

%% Space discretization

Data.degree        = [1, 2, 3, 4, 5, 6, 7, 8];		% Polynomial degree
Data.penalty_coeff = 10;   			                % Penalty coefficient
Data.epsilon       = 0;

Data.NonLinearSolver = 'Semi-implicit';     % Nonlinearity treatment: Semi-implicit or Picard iterations

%% Nonlinear solver settings

Data.NLS_tolerance  = 1e-10;
Data.NLS_max_it     = 1000;
Data.NLS_relax      = 1;

%% Visualization settings

Data.PlotExact          = true;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 1;
Data.NPtsVisualization   = 5;
Data.parallel = false;            % Parallel settings
Data.cores = 1;                   % Number of cores used in the parallel setting

%% Assembling settings
Data.quadrature = 'QF';           % Method of integrals calculation: Quadrature-Free - Subtriangulation

%% Choice of the formulation
Data.formulation = 'IPDG';           % Discrete formulation choice: IPDG - LDG - SPLDG

%% Adaptivity
Data.Adaptivity         = false;                % Flag of adaptivity                                             
Data.maxDegree          = 2;                    % Maximum polynomial degree
Data.AdaptFunc          = @(tau_r) floor(1 + 2*Data.maxDegree/pi* atan(tau_r));     % Adaptivity function
Data.AdaptIts           = 1;                   % Maximum iterations of adaptivity
Data.AdaptivityStep     = 1; 
