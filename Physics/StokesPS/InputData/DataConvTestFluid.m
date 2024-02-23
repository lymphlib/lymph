%% Unsteady Stokes problem in [0 1]^2 with mixed conditions

Data.name = 'DataTestFluid';

Data.TagElFluid   = 1;         % Element tag
Data.TagBcFluid   = [2 3 4 5]; % Boundary tag
Data.LabBcFluid   = 'NDDN';    % (D)irichlet (N)eumann 
% REMEMBER sigma.n    = gn is Dirichlet
%          div(sigma) = gd is Neumann 

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 50;        % number of elements for a new mesh
Data.MeshFromFile = true;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = ["UnitSquareMixed_50_el", "UnitSquareMixed_100_el", ...
                     "UnitSquareMixed_200_el", "UnitSquareMixed_400_el"]; %filename for mesh 


%% Discretization properties                            
%% Time integration
Data.t0 = 0;
Data.T  =  1;
Data.dt = 0.001;

Data.timeint   = 'CN';

%% Space discretization
Data.degree  = 1;   % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Visualization settings
Data.PlotExact   = true;
Data.PlotGridSol = true;
Data.VisualizationStep = 0.1;
Data.PlotIniCond  = false;
Data.ComputeVelAndPres = false;

%% Save Solution settings
Data.timesave  = 1000;

%% properties fluid material
% Material parameter
Data.mu_f  = {@(x,y) 1 + 0.*x.*y};

% Forcing term fluid media sigma_11, sigma_12, sigma_21, sigma_22
Data.source_sigma   = {@(x,y) (-3 + 0.*x.*y); 
                       @(x,y)  0.*x.*y;
                       @(x,y)  0.*x.*y;
                       @(x,y)  0.*x.*y};
Data.source_sigma_t = {@(t) sin(2*t)};
% Data.source_sigma_t = {@(t) 1+0*t};

Data.source_sigma_d   = {@(x,y)  x.^2/2 + y.^2/4; 
                         @(x,y)  x.*y;
                         @(x,y)  x.*y;
                         @(x,y)  -x.^2/2 - y.^2/4};
Data.source_sigma_d_t = {@(t) 2*cos(2*t)};
% Data.source_sigma_d_t = {@(t) 0.*t};

% Dirichlet BC
Data.DirBCsigma      = {@(x,y)  x.^2;
                        @(x,y)  x.*y;
                        @(x,y)  x.*y;
                        @(x,y)  - 0.5*y.^2};                     
% Neumann BC
Data.NeuBCsigma_x  = {@(x,y)   2*x; 
                      @(x,y)   y;
                      @(x,y)   y;
                      @(x,y)   0.*x.*y};
 
Data.NeuBCsigma_y  = {@(x,y)   0.*x.*y; 
                      @(x,y)   x;
                      @(x,y)   x;
                      @(x,y)   -y};

                  
% exact solution --> used to compute the initial conditions
Data.sigma_ex    = {@(x,y)   x.^2;
                    @(x,y)   x.*y;
                    @(x,y)   x.*y;
                    @(x,y)   - 0.5*y.^2};
Data.sigma_t_ex  = {@(t) sin(2*t)}; 
Data.sigma_dt_ex = {@(t) 2*cos(2*t)};
% Data.sigma_t_ex  = {@(t) 1}; 
% Data.sigma_dt_ex = {@(t) 0.*t};

% forcing term for velocity and pressure recovery
Data.source_vel   = {@(x,y) 0.*x.*y; 
                     @(x,y) 0.*x.*y};
Data.source_vel_t = {@(t) 0*t};

Data.source_vel_d   = {@(x,y) 0.*x.*y; 
                       @(x,y) 0.*x.*y};
Data.source_vel_d_t = {@(t) 0*t};

Data.vel0 = {@(x,y)     0.*x.*y; 
             @(x,y)     0.*x.*y};


