%% Elastodynamics problem in [0 1]^2 with Dirichlet conditions

Data.name = 'DataTestEla';

Data.TagElEla   = 1;         % Element tag
Data.TagBcEla   = [2 3 4 5]; % Boundary tag
Data.LabBcEla   = 'DDDD';    % (D)irichlet/(N)eumann/(A)bso

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 100;        % number of elements for a new mesh
Data.MeshFromFile = false;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'UnitSquare_200_el'; %filename for mesh 

%% Discretization properties 

% Time integration
Data.t0      = 0;
Data.T       = 0.01;
Data.dt      = 0.01;
Data.timeint = 'newmark';
Data.BetaNM  = 0.25;
Data.GammaNM = 0.5;

% Space discretization
Data.degree        = 2;  % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Quadrature settings
Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = true;
Data.PlotGridSol       = false;
Data.VisualizationStep = 1;
Data.PlotIniCond       = false;
Data.NqnVisualization  = 3; % Data.NqnVisualization must be odd and strictly greater than 1 (see Gauleg.m) 

%% Properties elastic material
Data.rho_el    = {@(x,y) 1 + 0.*x.*y};
Data.vs_el     = {@(x,y) 1 + 0.*x.*y};
Data.vp_el     = {@(x,y) 2 + 0.*x.*y};
Data.zeta      = {@(x,y) 0 + 0.*x.*y};

Data.mu_el     = {@(x,y) 1 + 0.*x.*y}; % Data.vs_el^2 * Data.rho_el;
Data.lam_el    = {@(x,y) 2 + 0.*x.*y}; % Data.vp_el^2 * Data.rho_el - 2*Data.mu_el;

%% Forcing terms
% forcing term elastic media
Data.source_f   = {@(x,y,t) sin(sqrt(2)*pi*t).*(2*pi^2*cos(pi*y).*sin(pi*y).*(4*cos(2*pi*x) - (2*cos(2*pi*x))/2 - 1)); ... % u_e
                   @(x,y,t) sin(sqrt(2)*pi*t).*(-2*pi^2*cos(pi*x).*sin(pi*x).*(4*cos(2*pi*y) - (2*cos(2*pi*y))/2 - 1))};

Data.source_g   = {@(x,y,t) 0*t.*x.*y; ... % \dot{u}_e
                   @(x,y,t) 0*t.*x.*y};

% Moment source tensor
Data.sourceMxx_el = {@(x,y,t) 0*t.*x.*y}; % Forcing term Mxx
Data.sourceMyy_el = {@(x,y,t) 0*t.*x.*y}; % Forcing term Myy
Data.sourceMxy_el = {@(x,y,t) 0*t.*x.*y}; % Forcing term Mxy
Data.sourceMyx_el = {@(x,y,t) 0*t.*x.*y}; % Forcing term Myx

%% Boundary conditions
% Dirichlet boundary conditions
Data.DirBCEla = {@(x,y,t) sin(sqrt(2)*pi*t).*(-sin(pi*x).^2 .* sin(2*pi*y)); ...
                 @(x,y,t) sin(sqrt(2)*pi*t).*( sin(pi*y).^2 .* sin(2*pi*x))};

%% Initial conditions (& exact solutions for errors computation)
% Exact solution -> used to compute the initial conditions
Data.ue_ex = {@(x,y,t) sin(sqrt(2)*pi*t).*(-sin(pi*x).^2 .* sin(2*pi*y)); ...
              @(x,y,t) sin(sqrt(2)*pi*t).*( sin(pi*y).^2 .* sin(2*pi*x))};

% Time-derivative of the exact solution -> used to compute the initial conditions
Data.due_dt_ex = {@(x,y,t) sqrt(2)*pi*cos(sqrt(2)*pi*t).*(-sin(pi*x).^2 .* sin(2*pi*y)); ...
                  @(x,y,t) sqrt(2)*pi*cos(sqrt(2)*pi*t).*( sin(pi*y).^2 .* sin(2*pi*x))};

% Exact gradient -> used for the error analysis
Data.grad_ue_ex =  {@(x,y,t) sin(sqrt(2)*pi*t).*(-2*pi*cos(pi*x).*sin(pi*x).*sin(2*pi*y)); ... % du1/dx
                    @(x,y,t) sin(sqrt(2)*pi*t).*(-2*pi*cos(2*pi*y).*sin(pi*x).^2); ...         % du1/dy
                    @(x,y,t) sin(sqrt(2)*pi*t).*( 2*pi*cos(2*pi*x).*sin(pi*y).^2); ...         % du2/dx
                    @(x,y,t) sin(sqrt(2)*pi*t).*( 2*pi*cos(pi*y).*sin(2*pi*x).*sin(pi*y))};    % du2/dy
