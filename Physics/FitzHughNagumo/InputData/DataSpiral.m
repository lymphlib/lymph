%% Simulation for FHN equation - Spiral wave

Data.name = 'FHN';

Data.LabEl = {'N'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'NNNN'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 0;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain          = [0 2.5 0 2.5];
Data.N               = [50 50]; 
Data.meshfromfile    = false;
Data.foldername      = 'InputMesh';
Data.VTKMeshFileName = 'WavesMesh.vtk';
Data.meshfileseq     = "NeuFHN_2500_el.mat";

Data.isotropy    = true;
Data.AxnDiffFile = "...";

%% Material properties 
d     = 1e-4;

Data.D_ext       = {@(x,y,t) d+0.*x};  % External diffusion
Data.D_axn       = {@(x,y,t) 0.*x};       % Axonal diffusion
Data.Cm          = {@(x,y,t) 1 + 0.*x}; 
Data.Chi         = {@(x,y,t) 1 + 0.*x};

Data.k       = {@(x,y,t) 1 + 0.*x};
Data.epsilon = {@(x,y,t) 0.01 + 0.*x};
Data.gamma   = {@(x,y,t) 1 + 0.*x};
Data.beta    = {@(x,y,t) 0.5 + 0.*x};
Data.a       = {@(x,y,t) 0.1 + 0.*x};

% Forcing Term
Data.homog_source_f = true;
Data.source_f   = {@(x,y,t) 0.*x};

% Boundary Conditions
Data.DirBC    =  {@(x,y,t)  0.*x};
Data.gradDirBC =  {@(x,y,t) 0.*x};

% Exact Solution (if any)
Data.u_ex = {@(x,y,t) double( ...
    (x <= 1.25) & ...
    (y <= 1.25) )};

eps = 0.02; 
Data.w_ex = {@(x,y,t) 0.1 * 0.5 * (1 + tanh((y - 1.25)/eps))};

% Gradient of the Exact Solution
Data.du_dx_ex = {@(x,y,t) 0.*pi*cos(pi*x).*sin(pi*y)*exp(-t)};
Data.du_dy_ex = {@(x,y,t) 0.*pi*sin(pi*x).*cos(pi*y)*exp(-t)};
Data.du_dt_ex =  {@(x,y,t) -0.*sin(pi*x).*sin(pi*y)*exp(-t)};

%% Discretization properties

%% Time discretization
Data.t0    = 0;
Data.T     = 5000.00;
Data.dt    = 5e-1;

Data.theta    = 0.5;

%% Space discretization
Data.degree        = 4;		        % Polynomial degree
Data.penalty_coeff = 10;   			% Penalty coefficient

%% Nonlinear solver settings
Data.NonLinearSolver = 'Picard iterations';     % Nonlinearity treatment: Semi-implicit or Picard iterations

Data.NLS_tolerance  = 1e-10;
Data.NLS_max_it     = 1000;
Data.NLS_relax      = 1;

%% Visualization settings
Data.PlotExact          = false;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 40;
Data.NPtsVisualization  = 5;

%% Assembling settings
Data.quadrature = 'QF';           % Method of integrals calculation: Quadrature-Free - Subtriangulation

%% Adaptivity
Data.Adaptivity         = true;                % Flag of adaptivity                                             
Data.maxDegree          = 4;                    % Maximum polynomial degree
Data.AdaptFunc          = @(tau_r) floor(1 + (Data.maxDegree-1)*tanh(5*tau_r));     % Adaptivity function
Data.AdaptIts           = 2;                   % Maximum iterations of adaptivity
Data.AdaptivityStep     = 10; 
