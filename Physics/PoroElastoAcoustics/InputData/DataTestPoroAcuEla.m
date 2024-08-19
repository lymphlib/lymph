% Domain Poro --> (-1,0)x(0,1) U (0,1)x(0,1) <-- Ela 
Data.name = 'DataConvTestPoroAcuEla';
Data.NPhys = 3;

Data.TagElPoro  = 1; % Element tag
Data.TagBcPoro  = [12 5]; % Boundary tag
Data.LabBcPoro  = 'DD'; % (D)irichlet/(N)eumann/(A)bso

Data.TagElAcu   = [2 3]; % Element tag
Data.TagBcAcu   = [8 9 10 11]; % Boundary tag
Data.LabBcAcu   = 'DDDD'; % (D)irichlet/(N)eumann 

Data.TagElEla   = 4; % Element tag
Data.TagBcEla   = [6 7]; % Boundary tag
Data.LabBcEla   = 'DD'; % (D)irichlet/(N)eumann/(A)bso

%% Geometrical properties 
Data.domain       = [-1 1 -1 1]; % domain bounds for a new mesh
Data.N            = 2000;        % number of elements for a new mesh
Data.MeshFromFile = false;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'Quadridomain'; %filename for mesh 

%% Discretization properties                            
%% Time integration
Data.t0 = 0;
Data.T  =  0.001;
Data.dt = 0.001;

Data.timeint   = 'newmark';
Data.BetaNM = 0.25;
Data.GammaNM = 0.5;

%% Space discretization
Data.degree  = 3;   % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Quadrature settings
Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = false;
Data.PlotGridSol       = true;
Data.VisualizationStep = 0.1;
Data.PlotIniCond       = false;
Data.NqnVisualization = 3;

%% Save Solution settings
Data.VisualizationStep  = 1100;


%% properties porous-materials
Data.rho_f   = {@(x,y) 1 + 0.*x.*y};
Data.rho_s   = {@(x,y) 1 + 0.*x.*y};
Data.phi_por = {@(x,y) 0.5 + 0.*x.*y};
Data.zetap   = {@(x,y) 0.*x.*y};
Data.mu      = {@(x,y) 1 + 0.*x.*y};
Data.lam     = {@(x,y) 1 + 0.*x.*y};
Data.beta    = {@(x,y) 1 + 0.*x.*y};
Data.m       = {@(x,y) 1 + 0.*x.*y};
Data.eta     = {@(x,y) 0 + 0.*x.*y};
Data.k_per   = {@(x,y) 1 + 0.*x.*y};
Data.a_coef  = {@(x,y) 1 + 0.*x.*y};
Data.tau     = 0; % fluid flow poro-acoustic
Data.delta   = 0; % fluid flow poro-elastic

% ATT !
% Data.rho     = Data.phi_por .* Data.rho_f + ( 1-Data.phi_por ) .* Data.rho_s;
% Data.rho_w   = (Data.a_coef ./ Data.phi_por) .* Data.rho_f;
% Data.vs_poro = sqrt(Data.mu .* Data.rho_w./(Data.rho.*Data.rho_w-Data.rho_f.^2));

Data.rho     = {@(x,y) 1 + 0.*x.*y};
Data.rho_w   = {@(x,y) 2 + 0.*x.*y};
Data.vs_poro = {@(x,y) sqrt(2) + 0.*x.*y};

% for i = 1 : length(Data.rho)
%     A = [Data.rho(i) Data.rho_f(i); Data.rho_f(i) Data.rho_w(i)];
%     B = [Data.lam(i) + 2*Data.mu(i) + Data.m(i)*Data.beta(i)^2,  Data.m(i)*Data.beta(i);
%         Data.m(i)*Data.beta(i),                                  Data.m(i)];
%     L = sqrt(eig(A^-1*B));
%     Data.vp_poroI(i)  = max(L);
%     Data.vp_poroII(i) = min(L);
% end

Data.vp_poroI = {@(x,y)  2.557612414958355 + 0.*x.*y};
Data.vp_poroII = {@(x,y) 0.677213950573148 + 0.*x.*y};


% forcing term porous media u_p
Data.source_up = {@(x,y) -((8*x.^2.*cos((pi*x)/2).^2.*sin((pi*x)/2) + 24*y.^2.*cos((pi*x)/2).^2.*sin((pi*x)/2) + 6*x.^2.*y.^2*pi^2.*sin((pi*x)/2) - 48*x.*y.^2*pi.*cos((pi*x)/2) - 16*x.^2.*y*pi.*cos((pi*x)/2) + 32*x.*y.*cos((pi*x)/2).^2.*sin((pi*x)/2) + 72*x.*y.^2*pi.*cos((pi*x)/2).^3 + 24*x.^2.*y*pi.*cos((pi*x)/2).^3 - 27*x.^2.*y.^2*pi^2.*cos((pi*x)/2).^2.*sin((pi*x)/2)))/2; ...
                  @(x,y) -((24*x.^2.*cos((pi*x)/2).^2.*sin((pi*x)/2) + 8*y.^2.*cos((pi*x)/2).^2.*sin((pi*x)/2) + 2*x.^2.*y.^2*pi^2.*sin((pi*x)/2) - 16*x.*y.^2*pi.*cos((pi*x)/2) - 16*x.^2.*y*pi.*cos((pi*x)/2) + 32*x.*y.*cos((pi*x)/2).^2.*sin((pi*x)/2) + 24*x.*y.^2*pi.*cos((pi*x)/2).^3 + 24*x.^2.*y*pi.*cos((pi*x)/2).^3 - 9*x.^2.*y.^2*pi^2.*cos((pi*x)/2).^2.*sin((pi*x)/2)))/2};
Data.source_up_t = {@(t) cos(pi*t*sqrt(2))};

% forcing term porous media \dot{u}_p
Data.source_upd = {@(x,y) 0.*x.*y; ...
                   @(x,y) 0.*x.*y};
Data.source_upd_t = {@(t) 0.*t};

% forcing term porous media w_p
Data.source_wp = {@(x,y) 2*x.^2.*y.^2*pi^2.*cos((pi*x)/2).*sin(pi*x); ...
                  @(x,y) 2*x.^2.*y.^2*pi^2.*cos((pi*x)/2).*sin(pi*x)};
Data.source_wp_t = {@(t) cos(pi*t*sqrt(2))};

% forcing term porous media \dot{w}_p
Data.source_wpd = {@(x,y) 0.*x.*y; ...
                   @(x,y) 0.*x.*y};
Data.source_wpd_t = {@(t) 0.*t};

% Moment source tensor
Data.sourceMxx_poro   = {@(x,y) 0*x.*y}; % Forcing term Mxx
Data.sourceMyy_poro   = {@(x,y) 0*x.*y}; % Forcing term Myy
Data.sourceMxy_poro   = {@(x,y) 0*x.*y}; % Forcing term Mxy
Data.sourceMyx_poro   = {@(x,y) 0*x.*y}; % Forcing term Myx

% Dirichlet BC
Data.DirBCPoro_up   = {@(x,y) y.^2.*x.^2.*cos(pi*x/2).*sin(pi*x); ...
                       @(x,y) y.^2.*x.^2.*cos(pi*x/2).*sin(pi*x)};
Data.DirBCPoro_up_t = {@(t) cos(pi*t*sqrt(2))};

Data.DirBCPoro_wp = {@(x,y) -y.^2.*(x.^2.*cos(pi*x/2).*sin(pi*x)); ...
                     @(x,y) -y.^2.*(x.^2.*cos(pi*x/2).*sin(pi*x))};
Data.DirBCPoro_wp_t = {@(t) cos(pi*t*sqrt(2))};

% exact solution --> used to compute the initial conditions and Dirichlet
% (time dependence)
Data.up_ex    =  {@(x,y) y.^2.*x.^2.*cos(pi*x/2).*sin(pi*x); ...
                  @(x,y) y.^2.*x.^2.*cos(pi*x/2).*sin(pi*x)};
Data.up_t_ex  =  {@(t) cos(pi*t*sqrt(2))};
Data.dup_t_ex =  {@(t) -sqrt(2)*pi*sin(sqrt(2)*pi*t)};

Data.wp_ex    =  {@(x,y) -y.^2.*(x.^2.*cos(pi*x/2).*sin(pi*x)); ...
                  @(x,y) -y.^2.*(x.^2.*cos(pi*x/2).*sin(pi*x))};
Data.wp_t_ex  =  {@(t) cos(pi*t*sqrt(2))};
Data.dwp_t_ex =  {@(t) -sqrt(2)*pi*sin(sqrt(2)*pi*t)};

% exact gradient --> used for the error analysis
% du1/dx, du1/dy, du2/dx, du2/dy
% Data.grad_up_ex =  {@(x,y) 2*x.*cos((pi*x)/2).*sin(pi*x) + x.^2*pi.*cos(pi*x).*cos((pi*x)/2) - (x.^2*pi.*sin(pi*x).*sin((pi*x)/2))/2; ...
%                    @(x,y) 0*x.*y; ...
%                    @(x,y) 2*x.*cos((pi*x)/2).*sin(pi*x) + x.^2*pi*cos(pi*x).*cos((pi*x)/2) - (x.^2*pi*sin(pi*x).*sin((pi*x)/2))/2; ...
%                    @(x,y) 0*x.*y;};

%Data.grad_wp_ex =  {@(x,y) -(2.*x*cos((pi*x)/2).*sin(pi*x) + x.^2*pi*cos(pi*x).*cos((pi*x)/2) - (x.^2*pi*sin(pi*x).*sin((pi*x)/2))/2); ...
%                    @(x,y) 0*x.*y; ...
%                    @(x,y) -(2.*x.*cos((pi*x)/2).*sin(pi*x) + x.^2*pi*cos(pi*x).*cos((pi*x)/2) - (x.^2*pi*sin(pi*x).*sin((pi*x)/2))/2); ...
%                    @(x,y) 0*x.*y;};

%% properties acoustic material

Data.rho_a = {@(x,y) 1 + 0.*x.*y; @(x,y) 1 + 0.*x.*y};
Data.c     = {@(x,y) 1 + 0.*x.*y; @(x,y) 1 + 0.*x.*y};

% forcing term acoustic media
Data.source_phi   = {@(x,y) -(2*x.^2.*sin(pi*x).*sin(pi*y) + 2*y.^2.*sin(pi*x).*sin(pi*y) + 4*x.*y.^2*pi.*cos(pi*x).*sin(pi*y) + 4*x.^2.*y*pi.*cos(pi*y).*sin(pi*x) - 2*x.^2.*y.^2*pi^2.*sin(pi*x).*sin(pi*y) + 2*x.^2.*y.^2*pi^2.*sin(pi*x).*sin(pi*y))};
Data.source_phi_t = {@(t) sin(pi*t*sqrt(2))};

% Dirichlet conditions
Data.DirBCAcu   = {@(x,y) y.^2.*x.^2.*sin(pi*x).*sin(pi*y)};
Data.DirBCAcu_t = {@(t) sin(pi*t*sqrt(2))};

% exact solution --> used to compute the initial conditions
Data.phi_ex    = {@(x,y) y.^2.*x.^2.*sin(pi*x).*sin(pi*y)};
Data.phi_t_ex  = {@(t) sin(pi*t*sqrt(2))}; 
Data.dphi_t_ex = {@(t) sqrt(2)*pi*cos(sqrt(2)*pi*t)};

% exact gradient --> used for the error analysis
% Data.grad_phi_ex =  {@(x,y) 2*x.*sin(pi*x).*sin(pi*y) + x.^2*pi*cos(pi*x).*sin(pi*y); ...
%                      @(x,y) x.^2*pi*cos(pi*y).*sin(pi*x)};



%% Properties elastic material
Data.rho_el    = {@(x,y) 1 + 0.*x.*y};
Data.vs_el     = {@(x,y) 1 + 0.*x.*y};
Data.vp_el     = {@(x,y) 2 + 0.*x.*y};
Data.zeta      = {@(x,y) 0 + 0.*x.*y};

Data.mu_el     = {@(x,y) 1 + 0.*x.*y}; % Data.vs_el^2 * Data.rho_el;
Data.lam_el    = {@(x,y) 2 + 0.*x.*y}; % Data.vp_el^2 * Data.rho_el - 2*Data.mu_el;

% forcing term elastic media
Data.source_ue   = {@(x,y) -2*(x.^2.*sin(2*pi*x) + 4*y.^2.*sin(2*pi*x) + 12.*x.*y.*cos(2*pi*x).*sin(2*pi*x) + 16*x.*y.^2*pi.*cos(2*pi*x) + 12*x.^2.*y*pi.*(2*cos(2*pi*x).^2 - 1)); ...
                    @(x,y) -2*(4*x.^2.*sin(4*pi*x) + y.^2.*sin(4*pi*x) + 6.*x.*y.*sin(2*pi*x) + 6*x.^2.*y*pi.*cos(2*pi*x) + 8*x.*y.^2*pi.*cos(4*pi*x))};
Data.source_ue_t = {@(t) cos(4*pi*t)};

% forcing term elastic media \dot{u}_e
Data.source_ued   = {@(x,y) 0*x.*y; @(x,y) 0*x.*y};
Data.source_ued_t = {@(t) 0.*t};

% Moment source tensor
Data.sourceMxx_el   = {@(x,y) 0*x.*y}; % Forcing term Mxx
Data.sourceMyy_el   = {@(x,y) 0*x.*y}; % Forcing term Myy
Data.sourceMxy_el   = {@(x,y) 0*x.*y}; % Forcing term Mxy
Data.sourceMyx_el   = {@(x,y) 0*x.*y}; % Forcing term Myx

% Dirichlet BC
Data.DirBCEla  = {@(x,y) x.^2.*y.^2.*(sin(2*pi*x)); ...
                  @(x,y) x.^2.*y.^2.*(sin(4*pi*x))};

Data.DirBCEla_t = {@(t) cos(4*pi*t)};

% exact solution --> used to compute the initial conditions and Dirichlet
% (time dependence)
Data.ue_ex    =  {@(x,y) x.^2.*y.^2.*(sin(2*pi*x)); ...
                  @(x,y) x.^2.*y.^2.*(sin(4*pi*x))};
Data.ue_t_ex  =  {@(t) cos(4*pi*t)};
Data.due_t_ex =  {@(t) -4*pi*sin(4*pi*t)};

% exact gradient --> used for the error analysis
% du1/dx, du1/dy, du2/dx, du2/dy
% Data.grad_ue_ex =  {@(x,y) 2*x.*sin(2*pi*x) + 2*x.^2*pi.*cos(2*pi*x); ...
%                     @(x,y) 0.*x.*y; ...
%                     @(x,y) 2*x.*sin(2*pi*x) + 2*x.^2*pi.*cos(2*pi*x); ...
%                     @(x,y) 0.*x.*y};
