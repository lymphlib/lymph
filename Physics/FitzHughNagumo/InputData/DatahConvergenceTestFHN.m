%% Convergence test for FHN equation - h convergence 

Data.name = 'ConvTestDir';

Data.LabEl = {'N'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'DDDD'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 1;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain = [-1 1 -1 1];
Data.N      = 30;
Data.meshfromfile = true;
Data.foldername = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq = ["FHN_180_el.mat","FHN_300_el.mat", ... 
                     "FHN_800_el.mat","FHN_1500_el.mat"];

Data.isotropy = true;
Data.AxnDiffFile = "...";

%% Material properties 

Data.D_ext       = {@(x,y,t) 1+0.*x};     % External diffusion
Data.D_axn       = {@(x,y,t) 0+0.*x};     % Axonal diffusion

Data.Cm          = {@(x,y,t) 1 + 0.*x}; 
Data.Chi         = {@(x,y,t) 1 + 0.*x};

Data.k       = {@(x,y,t) 1 + 0.*x};
Data.epsilon = {@(x,y,t) 1 + 0.*x};
Data.gamma   = {@(x,y,t) 2 + 0.*x};
Data.beta    = {@(x,y,t) 1 + 0.*x};
Data.a       =  {@(x,y,t) 0.5 + 0.*x};

% Forcing Term
Data.homog_source_f = false;
Data.source_f =  {@(x,y,t) - sin(pi*x).*sin(pi*y)*exp(-t) + 2.*pi.*pi.*sin(pi*x).*sin(pi*y)*exp(-t) + ...
                            (sin(pi*x).*sin(pi*y)*exp(-t) - 1).*(sin(pi*x).*sin(pi*y)*exp(-t)- 0.5).*sin(pi*x).*sin(pi*y).*exp(-t) + ...
                            sin(pi*x).*sin(pi*y).*exp(-t)};
% Boundary Conditions
Data.DirBC    =  {@(x,y,t)  sin(pi*x).*sin(pi*y)*exp(-t)};
Data.gradDirBC =  {@(x,y,t) -pi*sin(pi*x).*cos(pi*y).*exp(-t), @(x,y,t) -pi*cos(pi*x).*sin(pi*y).*exp(-t)};

% Exact Solution (if any)
Data.u_ex     =  {@(x,y,t)  sin(pi*x).*sin(pi*y)*exp(-t)};
Data.w_ex     =  {@(x,y,t)  sin(pi*x).*sin(pi*y)*exp(-t)};

% Gradient of the Exact Solution
Data.du_dx_ex = {@(x,y,t) pi*cos(pi*x).*sin(pi*y)*exp(-t)};
Data.du_dy_ex = {@(x,y,t) pi*sin(pi*x).*cos(pi*y)*exp(-t)};
Data.du_dt_ex =  {@(x,y,t) -sin(pi*x).*sin(pi*y)*exp(-t)};


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
