%% Unsteady Stokes problem: flow around a circular cylinder

Data.name = 'DataTestFACFluid';

Data.TagElFluid   = 1; % Element tag
Data.TagBcFluid   = [2 3 4 5 6]; % Boundary tag
Data.LabBcFluid   = 'DNDDD'; % Dirichlet/Neumann/Abso 
% REMEMBER sigma.n = gn is Dirichlet
%          div(sigma) = gd is Neumann 

%% Geometrical properties 
Data.domain       = [-1 4 -1 1]; % domain bounds for a new mesh
Data.N            = 50;        % number of elements for a new mesh
Data.MeshFromFile = true;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh_FAC.vtk';
Data.meshfileseq  = 'FlowAroundCilinder_2000_el'; %filename for mesh 

%% Discretization properties                            
%% Time integration
Data.t0 = 0;
Data.T  =  1;
Data.dt = 0.01;
Data.timeint   = 'CN';

%% Space discretization
Data.degree  = 3;   % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Visualization settings
Data.PlotExact   = false;
Data.PlotGridSol = false;
Data.VisualizationStep = 0.1;
Data.PlotIniCond  = false;
Data.ComputeVelAndPres = true;

%% Save Solution settings
Data.timesave  = 100;

%% properties fluid material
% Material parameter
Data.mu_f  = {@(x,y) 2 + 0.*x.*y};

% forcing term fluid media
% s_11, s_12, s_21, s22
Data.source_sigma   = {@(x,y)  0.*x.*y; 
                       @(x,y)  0.*x.*y;
                       @(x,y)  0.*x.*y;
                       @(x,y) 0.*x.*y};
Data.source_sigma_t = {@(t) 1+0.*t};

Data.source_sigma_d   = {@(x,y)  0.*x.*y; 
                         @(x,y)  0.*x.*y;
                         @(x,y)  0.*x.*y;
                         @(x,y) 0.*x.*y};
Data.source_sigma_d_t = {@(t) 0*t};


% REMEMBER sigma.n = gn is Dirichlet
%          div(sigma) = gd is Neumann 
% Dirichlet BC
Data.DirBCsigma      = {@(x,y)  0.*x.*y;
                        @(x,y)  0.*x.*y;
                        @(x,y)  0.*x.*y;
                        @(x,y)  0.*x.*y};
                     
% Neumann BC
Data.NeuBCsigma_x  = {@(x,y)  (1-y.^2).*((x+1)<=1.e-3) + 0.*x; 
                      @(x,y)  0.*x.*y;
                      @(x,y)  0.*x.*y;
                      @(x,y)  0.*x.*y};
 
Data.NeuBCsigma_y  = {@(x,y)  0.*x.*y; 
                      @(x,y)  0.*x.*y;
                      @(x,y)  0.*x.*y;
                      @(x,y)  0.*x.*y};

                  
% exact solution --> used to compute the initial conditions
Data.sigma_ex    = {@(x,y)  0.*x.*y;
                    @(x,y)  0.*x.*y;
                    @(x,y)  0.*x.*y;
                    @(x,y)  0.*x.*y};
                     
Data.sigma_t_ex  = {@(t) 0*t}; 
Data.sigma_dt_ex = {@(t) 0*t};

% forcing term for velocity and pressure recovery
Data.source_vel   = {@(x,y)     0.*x.*y; 
                     @(x,y)     0.*x.*y};
Data.source_vel_t = {@(t) 0.*t};

Data.source_vel_d   = {@(x,y)  0.*x.*y ; 
                       @(x,y)  0.*x.*y};
Data.source_vel_d_t = {@(t) 0*t};

Data.vel0 = {@(x,y)     0.*x.*y; 
             @(x,y)     0.*x.*y};
