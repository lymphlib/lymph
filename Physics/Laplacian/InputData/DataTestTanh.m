
%% Adaptive Poisson problem in [0,1]^2 with Dirichlet conditions

Data.name = 'AdaptLap';

Data.LabEl = {'L'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'DDDD'};  % (D)irichlet/(N)eumann/(A)bsorbing

%% Geometrical properties 
Data.domain       = [0 1 0 1]; % domain bounds for a new mesh
Data.N            = 180;     % number of elements for a new mesh
Data.MeshFromFile = false;      % read mesh from file
Data.FolderName   = 'InputMesh';
Data.VTKMeshFileName = 'Mesh.vtk';
Data.meshfileseq  = 'Lap'; %filename for mesh

%% Space discretization
Data.degree  = 6;        % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Quadrature settings
Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = true;
Data.PlotGridSol       = true;
Data.NPtsVisualization = 5;

%% Material properties 
Data.mu       = {@(x,y) 1.*x.^0.*y.^0};

% Forcing Term
Data.source   = {@(x,y) -80*(sech(20*(-0.8+x.^2+y.^2))).^2.*(-1+40*(x.^2+y.^2).*tanh(20*(-0.8+x.^2+y.^2)))};

% Boundary Conditions
Data.DirBC    = {@(x,y) tanh(-20*(x.^2+y.^2-0.8))};

% Gradient of the Exact Solution
Data.gradDirBC =  {@(x,y) -40*x.*sech(20.*(x.^2+y.^2-4/5)).^2; ...
                   @(x,y) -40*y.*sech(20.*(x.^2+y.^2-4/5)).^2};

% Exact Solution (if any)
Data.u_ex     = {@(x,y) tanh(-20*(x.^2+y.^2-0.8))};

% Gradient of the Exact Solution
Data.gradu_ex =  {@(x,y) -40*x.*sech(20.*(x.^2+y.^2-4/5)).^2; ...
                  @(x,y) -40*y.*sech(20.*(x.^2+y.^2-4/5)).^2};

%% Adaptivity
Data.Adaptivity = true;                 % Flag of adaptivity                                             
Data.maxDegree  = 6;                    % Maximum polynomial degree
Data.AdaptFunc  = @(tau_r) floor(1 + 2*Data.maxDegree/pi* atan(tau_r));     % Adaptivity function
Data.AdaptIts   = 10;                   % Maximum iterations of adaptivity
