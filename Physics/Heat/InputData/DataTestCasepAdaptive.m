%% Patch test for Heat equation

Data.name = 'TestCase';


Data.LabEl = {'L'};     % Element labels
Data.TagEl = { 1 };     % Element tags

Data.TagBc = {[2 3 4 5]};  % Boundary tags
Data.LabBc = {'DDDD'};  % (D)irichlet/(N)eumann/(A)bsorbing

Data.TagApplyBCs = 1;      % Skip the assembling of BCs if not necessary

%% Geometrical properties

Data.domain             = [-1 1 -1 1];
Data.N                  = 500;                     % Number of mesh elements
Data.MeshFromFile       = true;
Data.FolderName         = 'InputMesh';
Data.VTKMeshFileName    = 'Mesh.vtk';
Data.meshfileseq        = 'SquareCentered_500_el.mat'; % Names of mesh files

%% Material properties 

Data.mu          = {@(x,y) 1+0.*x};     % Diffusion parameter
Data.sigma       = {@(x,y) 0+0.*x};       % Reaction parameter
epsilon = 0.08;

% Forcing Term
Data.homog_source_f = false;
Data.source_f       = {@(x,y,t) exp((1 - exp((x - t.*y)./epsilon))./(1 - exp(-1./epsilon))).*(y.*exp((x - t.*y)./epsilon)./(epsilon.*(1 - exp(-1./epsilon))) + (1 + t.^2).*exp((x - t.*y)./epsilon)./(epsilon.^2.*(1 - exp(-1./epsilon))) - (1 + t.^2).*exp(2.*(x - t.*y)./epsilon)./(epsilon.^2.*(1 - exp(-1./epsilon)).^2))};

% Boundary Conditions
Data.DirBC     =  {@(x,y,t) exp((1-exp((x-t.*y)/epsilon))./(1-exp(-1/epsilon)))};
Data.gradDirBC =  {@(x,y,t) - (1./epsilon).*exp((x - t.*y)./epsilon)./(1 - exp(-1./epsilon)).*exp((1 - exp((x - t.*y)./epsilon))./(1 - exp(-1./epsilon))), ...
                   @(x,y,t) (t./epsilon).*exp((x - t.*y)./epsilon)./(1 - exp(-1./epsilon)).*exp((1 - exp((x - t.*y)./epsilon))./(1 - exp(-1./epsilon)))};

% Exact Solution (if any)
Data.u_ex     =  {@(x,y,t) exp((1-exp((x-t.*y)/epsilon))./(1-exp(-1/epsilon)))};

% Gradient of the Exact Solution
Data.du_dx_ex =  {@(x,y,t) - (1./epsilon).*exp((x - t.*y)./epsilon)./(1 - exp(-1./epsilon)).*exp((1 - exp((x - t.*y)./epsilon))./(1 - exp(-1./epsilon)))};
Data.du_dy_ex =  {@(x,y,t) (t/epsilon) .* exp((x - t.*y)/epsilon) ./ (1 - exp(-1/epsilon)) .* exp( (1 - exp((x - t.*y)/epsilon)) ./ (1 - exp(-1/epsilon))) };
Data.du_dt_ex =  {@(x,y,t) (y/epsilon) .* exp((x - t.*y)/epsilon)./ (1 - exp(-1/epsilon)) .* exp( (1 - exp((x - t.*y)/epsilon)) ./ (1 - exp(-1/epsilon)) )};

%% Discretization properties

%% Time discretization

Data.t0     = 0;
Data.T      = 1;
Data.dt     = 1e-2;
Data.theta  = 0.5;

%% Space discretization

Data.degree        = 5;          % Polynomial degree
Data.penalty_coeff = 10;         % Penalty coefficient

%% Quadrature settings

Data.quadrature = "ST";       % Quadrature type: ST/QF

%% Visualization settings

Data.PlotExact          = true;
Data.PlotIniCond        = true;
Data.PlotGridSol        = true;
Data.VisualizationStep  = 10;
Data.NPtsVisualization  = 5; 

%% Adaptivity
Data.Adaptivity         = true;                 % Flag of adaptivity                                             
Data.maxDegree          = 5;                    % Maximum polynomial degree
Data.AdaptFunc          = @(tau_r) floor(1 + 2*Data.maxDegree/pi* atan(tau_r));     % Adaptivity function
Data.AdaptIts           = 1;                   % Maximum iterations of adaptivity
Data.AdaptivityStep     = 10; 