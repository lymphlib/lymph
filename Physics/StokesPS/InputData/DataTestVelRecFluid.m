%% Unsteady Stokes problem in [0 1]^2 with mixed conditions

Data.name = 'DataTestVelRecFluid';

Data.TagElFluid   = 1; % Element tag
Data.TagBcFluid   = [2 3 4 5]; % Boundary tag
Data.LabBcFluid   = 'DNDD'; % Dirichlet/Neumann/Abso 
% REMEMBER sigma.n = gn is Dirichlet
%          div(sigma) = gd is Neumann 

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 100;        % number of elements for a new mesh
Data.MeshFromFile = true;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'UnitSquareMixed_100_el'; %filename for mesh 

%% Discretization properties                            
%% Time integration
Data.t0 = 0;
Data.T  =  1;
Data.dt = 0.01;

Data.timeint   = 'CN';

%% Space discretization
Data.degree  = 2;   % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Visualization settings
Data.PlotExact   = true;
Data.PlotGridSol = true;
Data.VisualizationStep = 0.1;
Data.PlotIniCond  = false;
Data.ComputeVelAndPres = true;

%% Save Solution settings
Data.timesave  = 100;


%% properties fluid material
% Material parameter
Data.mu_f  = {@(x,y) 1 + 0.*x.*y};

% forcing term fluid media
% s_11, s_12, s_21, s22
Data.source_sigma   = {@(x,y) 0.*x.*y; 
                       @(x,y) 0.*x.*y;
                       @(x,y) 0.*x.*y;
                       @(x,y) 0.*x.*y};
Data.source_sigma_t = {@(t) t.^2};

Data.source_sigma_d   = {@(x,y)  -y + 0.*x; 
                         @(x,y) 1-x + 0.*y;
                         @(x,y)    0.*x.*y;
                         @(x,y)   y + 0.*x};
Data.source_sigma_d_t = {@(t) 2*t};


% Dirichlet BC
Data.DirBCsigma      = {@(x,y) 1-y + 0.*x;
                        @(x,y) 1-x + 0.*y;
                        @(x,y)  0.*x.*y;
                        @(x,y) 1+y + 0.*x};
                     
% Neumann BC (dsigma/dx)
Data.NeuBCsigma_x  = {@(x,y)       0.*x.*y; 
                      @(x,y)  -1 + 0.*x.*y;
                      @(x,y)       0.*x.*y;
                      @(x,y)       0.*x.*y};
% (dsigma/dy) 
Data.NeuBCsigma_y  = {@(x,y)  -1 + 0.*x.*y; 
                      @(x,y)       0.*x.*y;
                      @(x,y)       0.*x.*y;
                      @(x,y)   1 + 0.*x.*y};

                  
% exact solution --> used to compute the initial conditions
Data.sigma_ex    = {@(x,y) 1-y + 0.*x;
                    @(x,y) 1-x + 0.*y;
                    @(x,y)       0.*x.*y;
                    @(x,y) 1+y + 0.*x};
                     
Data.sigma_t_ex  = {@(t) t.^2}; 
Data.sigma_dt_ex = {@(t) 2.*t};



% forcing term for velocity and pressure recovery
Data.source_vel   = {@(x,y)     0.*x.*y; 
                     @(x,y) 1 + 0.*x.*y};
Data.source_vel_t = {@(t) -t.^2};

Data.source_vel_d   = {@(x,y)  (1-x).*y ; 
                       @(x,y)  0.5.*y.^2};
Data.source_vel_d_t = {@(t) 2*t};

Data.vel0 = {@(x,y)     0.*x.*y; 
             @(x,y)     0.*x.*y};



