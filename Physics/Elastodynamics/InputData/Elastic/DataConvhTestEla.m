%% Convergence Test for Elastodynamics

Data.name = 'DataConvTestEla';

Data.TagElEla   = 1; % Element tag
Data.TagBcEla   = [2 3 4 5]; % Boundary tag
Data.LabBcEla   = 'DDDD'; % Dirichlet/Neumann/Abso/


%% Geometrical properties 

Data.domain          = [0 1 0 1];
Data.N               = {100, 200, 400, 800};  %number of elements
Data.MeshFromFile    = false;
Data.FolderName      = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq     = ["UnitSquare_100_el.mat", ...
                        "UnitSquare_200_el.mat","UnitSquare_400_el.mat","UnitSquare_800_el.mat"]; %filename for mesh 

%% Discretization properties                            

% Time integration
Data.t0      = 0;
Data.T       = 1e-2;
Data.dt      = 0.001;
Data.timeint = 'newmark';
Data.BetaNM  = 0.25;
Data.GammaNM = 0.5;

% Space discretization
Data.degree        = 3;  % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient


%% Quadrature settings
Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = false;
Data.PlotGridSol       = false;
Data.VisualizationStep = 1;
Data.PlotIniCond       = false;
Data.NPtsVisualization = 3;

%% Properties elastic material
Data.rho_el    = {@(x,y) 1 + 0.*x.*y};
Data.vs_el     = {@(x,y) 1 + 0.*x.*y};
Data.vp_el     = {@(x,y) 2 + 0.*x.*y};
Data.zeta      = {@(x,y) 0 + 0.*x.*y};

Data.mu_el     = {@(x,y) 1 + 0.*x.*y}; % Data.vs_el^2 * Data.rho_el;
Data.lam_el    = {@(x,y) 2 + 0.*x.*y}; % Data.vp_el^2 * Data.rho_el - 2*Data.mu_el;

% forcing term elastic media
Data.source_f   = {@(x,y,t)  sin(sqrt(2)*pi*t).*(2*pi^2*cos(pi*y).*sin(pi*y).*(4*cos(2*pi*x) - (2*cos(2*pi*x))/2 - 1)); ...
                   @(x,y,t) sin(sqrt(2)*pi*t).*(-2*pi^2*cos(pi*x).*sin(pi*x).*(4*cos(2*pi*y) - (2*cos(2*pi*y))/2 - 1))};

Data.source_g   = {@(x,y,t) 0*t.*x.*y;...
                   @(x,y,t) 0*t.*x.*y};

% Moment source tensor
Data.sourceMxx_el   = {@(x,y,t) 0*t.*x.*y}; % Forcing term Mxx
Data.sourceMyy_el   = {@(x,y,t) 0*t.*x.*y}; % Forcing term Myy
Data.sourceMxy_el   = {@(x,y,t) 0*t.*x.*y}; % Forcing term Mxy
Data.sourceMyx_el   = {@(x,y,t) 0*t.*x.*y}; % Forcing term Myx

% Dirichlet BC
Data.DirBCEla  = {@(x,y,t) sin(sqrt(2)*pi*t).*(-sin(pi*x).^2 .* sin(2*pi*y)); ...
                  @(x,y,t) sin(sqrt(2)*pi*t).*( sin(pi*y).^2 .* sin(2*pi*x))};

% exact solution --> used to compute the initial conditions
Data.ue_ex = {@(x,y,t) sin(sqrt(2)*pi*t).*(-sin(pi*x).^2 .* sin(2*pi*y)); ...
              @(x,y,t) sin(sqrt(2)*pi*t).*( sin(pi*y).^2 .* sin(2*pi*x))};

Data.due_dt_ex = {@(x,y,t) sqrt(2)*pi*cos(sqrt(2)*pi*t).*(-sin(pi*x).^2 .* sin(2*pi*y)); ...
                  @(x,y,t) sqrt(2)*pi*cos(sqrt(2)*pi*t).*( sin(pi*y).^2 .* sin(2*pi*x))};

Data.ue_t_ex  =  {@(t) sin(sqrt(2)*pi*t)};
Data.due_t_ex =  {@(t) sqrt(2)*pi*cos(sqrt(2)*pi*t)};

% exact gradient --> used for the error analysis
% du1/dx, du1/dy, du2/dx, du2/dy
Data.grad_ue_ex =  {@(x,y,t) sin(sqrt(2)*pi*t).*(-2*pi*cos(pi*x).sin(pi*x).*sin(2*pi*y)); ...
                    @(x,y,t) sin(sqrt(2)*pi*t).*(-2*pi*cos(2*pi*y).*sin(pi*x).^2); ...
                    @(x,y,t) sin(sqrt(2)*pi*t).*( 2*pi*cos(2*pi*x).*sin(pi*y).^2); ...
                    @(x,y,t) sin(sqrt(2)*pi*t).*( 2*pi*cos(pi*y).*sin(2*pi*x).*sin(pi*y))};