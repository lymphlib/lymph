%> @file  MatPoroST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date 8 August 2024
%> @brief Assembly of the matrices for the poroelastic problem \cite ABNM2021
%>
%==========================================================================
%> @section classMatPoroPEAST Class description
%==========================================================================
%> @brief Assembly of the matrices for the poroelastic problem \cite ABNM2021
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Poroelastic matrices (Mass, Stiffness, dG, projection matrix, etc)
%>
%==========================================================================

function [Matrices] = MatPoroST(Data, neighbor, femregion)
%% Quadrature values
[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices
max_nedges = max(neighbor.nedges+1);

ii_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
jj_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
ii_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
jj_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

%% Initialization of the volume matrices

El.V1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.V2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.V3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.V4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.M1_P_rhof_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.M1_P_rhow_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.M1_P_rho_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.M1_P_eta_kper_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.M1_P_rho2_zeta_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.M1_P_rho_zeta2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.MPrjP_1_loc   = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
                  
El.B1_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B2_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B3_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B4_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.B1_beta_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B2_beta_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B3_beta_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B4_beta_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.B1_beta2_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B2_beta2_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B3_beta2_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.B4_beta2_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.P1_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.P2_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

El.P1_beta_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
El.P2_beta_loc  = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

%% Initialization of the edge matrices

Edge.S1_P_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_P_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.IT1_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.IT2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.IT3_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.IT4_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.S1_B_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S2_B_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S3_B_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_B_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.S1_B_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S2_B_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S3_B_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_B_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.S1_B_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S2_B_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S3_B_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_B_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.BT1_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT3_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT4_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.BT1_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT2_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT3_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT4_beta_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

Edge.BT1_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT2_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT3_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT4_beta2_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel_p),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel_p));

%% Initialization of the absorbing boundary matrices

Edge.ABC_uu_1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_uu_2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_uu_3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_uu_4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.ABC_uw_1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_uw_2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_uw_3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_uw_4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.ABC_wu_1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_wu_2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_wu_3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_wu_4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.ABC_ww_1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_ww_2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_ww_3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.ABC_ww_4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

%% Initialization of the acoustic boundary matrices

Edge.D1_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.D2_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.D3_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.D4_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
 
%% Initialization of the elastic boundary matrices

Edge.S1_P_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_P_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.S1_B_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S2_B_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S3_B_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_B_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.S1_B_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S2_B_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S3_B_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_B_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.S1_B_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S2_B_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S3_B_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.S4_B_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.BT1_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT3_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT4_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.BT1_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT2_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT3_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT4_beta_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

Edge.BT1_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT2_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT3_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));
Edge.BT4_beta2_EP_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel_p),femregion.nbases,femregion.nbases,ones(1,femregion.nel_p));

%% Loop over the elements

% Visualization of computational progress
index_shift = 0;
id_elas = max([Data.TagElAcu, Data.TagElPoro]);
if(isempty(id_elas)); id_elas = 0; end
nel_sh = 0;

prog = 0;
fprintf(1,'\t Computation Progress: %3d%%\n',prog);

for ie = 1 : femregion.nel

    ie_sh = ie - nel_sh;

    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    % Id and tag selection for elements
    id_ie  = femregion.id(ie);
    tag_ie = femregion.tag(ie);

    % Check if the element is poroelastic
    if tag_ie == 'P'

        % Selection of the matrix positions associated to element ie
        index = (ie_sh-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
        index_element = index_shift + (1:1:femregion.nedges(ie))';
        index_shift   = index_element(end);
        ii_index{ie_sh} = repmat(index, 1, femregion.nbases);
        jj_index{ie_sh} = repmat(index', femregion.nbases, 1);

        % Extraction of neighbor element and their edges
        neigh_ie      = neighbor.neigh{ie};
        neigh_ie_unq  = unique(neighbor.neigh{ie});
        neighedges_ie = neighbor.neighedges{ie};

        % Extraction of element geometrical information
        coords_ie          = femregion.coords_element{ie};
        [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);

        % Creation of the subtriangulation of the element
        edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
        Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
        Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);
        
        for iTria = 1:size(Tria,1)
            
            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);
            
            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            
            % Scaled weights
            dx = det(BJ) * w_2D;
            
            % Evaluation of physical parameters
            par.mu        = Data.mu{id_ie}(xq,yq);
            par.lam       = Data.lam{id_ie}(xq,yq);
            par.beta      = Data.beta{id_ie}(xq,yq);
            par.m         = Data.m{id_ie}(xq,yq);
            par.rho_f     = Data.rho_f{id_ie}(xq,yq);
            par.rho_w     = Data.rho_w{id_ie}(xq,yq);
            par.rho       = Data.rho{id_ie}(xq,yq);
            par.eta_kper  = Data.eta{id_ie}(xq,yq)./Data.k_per{id_ie}(xq,yq);
            par.rho2_zeta = 2 * Data.rho{id_ie}(xq,yq) .* Data.zetap{id_ie}(xq,yq);
            par.rho_zeta2 = Data.rho{id_ie}(xq,yq) .* Data.zetap{id_ie}(xq,yq).^2;

            % Construction and evalutation on the quadrature points of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);
                    
            % Local matrix assembling
            El.V1_loc{ie_sh} = El.V1_loc{ie_sh} + (dx .* ((par.lam+2*par.mu) .* gradqx))' * gradqx + (dx .* (par.mu .* gradqy))' * gradqy;
            El.V2_loc{ie_sh} = El.V2_loc{ie_sh} + (dx .* (par.lam .* gradqx))' * gradqy + (dx .* (par.mu .* gradqy))' * gradqx;
            El.V3_loc{ie_sh} = El.V3_loc{ie_sh} + (dx .* (par.lam .* gradqy))' * gradqx + (dx .* (par.mu .* gradqx))' * gradqy;
            El.V4_loc{ie_sh} = El.V4_loc{ie_sh} + (dx .* ((par.lam+2*par.mu) .* gradqy))' * gradqy + (dx .* (par.mu .* gradqx))' * gradqx;

            El.M1_P_rhof_loc{ie_sh} = El.M1_P_rhof_loc{ie_sh}  + (dx .* (par.rho_f .* phiq))' * phiq;
            El.M1_P_rhow_loc{ie_sh} = El.M1_P_rhow_loc{ie_sh}  + (dx .* (par.rho_w .* phiq))' * phiq;
            El.M1_P_rho_loc{ie_sh}  = El.M1_P_rho_loc{ie_sh}  + (dx .* (par.rho .* phiq))' * phiq;

            El.MPrjP_1_loc{ie_sh}   = El.MPrjP_1_loc{ie_sh}   + (dx .* (1         .* phiq))' * phiq;

            El.M1_P_rho2_zeta_loc{ie_sh} = El.M1_P_rho2_zeta_loc{ie_sh}  + (dx .* (par.rho2_zeta .* phiq))' * phiq;
            El.M1_P_rho_zeta2_loc{ie_sh} = El.M1_P_rho_zeta2_loc{ie_sh}  + (dx .* (par.rho_zeta2 .* phiq))' * phiq;
            El.M1_P_eta_kper_loc{ie_sh}  = El.M1_P_eta_kper_loc{ie_sh}  +  (dx .* (par.eta_kper  .* phiq))' * phiq;

            El.B1_loc{ie_sh} = El.B1_loc{ie_sh} + (dx .* (par.m .* gradqx))' * gradqx;
            El.B2_loc{ie_sh} = El.B2_loc{ie_sh} + (dx .* (par.m .* gradqx))' * gradqy;
            El.B3_loc{ie_sh} = El.B3_loc{ie_sh} + (dx .* (par.m .* gradqy))' * gradqx;
            El.B4_loc{ie_sh} = El.B4_loc{ie_sh} + (dx .* (par.m .* gradqy))' * gradqy;

            El.B1_beta_loc{ie_sh} = El.B1_beta_loc{ie_sh} + (dx .* (par.m .* par.beta .* gradqx))' * gradqx;
            El.B2_beta_loc{ie_sh} = El.B2_beta_loc{ie_sh} + (dx .* (par.m .* par.beta .* gradqx))' * gradqy;
            El.B3_beta_loc{ie_sh} = El.B3_beta_loc{ie_sh} + (dx .* (par.m .* par.beta .* gradqy))' * gradqx;
            El.B4_beta_loc{ie_sh} = El.B4_beta_loc{ie_sh} + (dx .* (par.m .* par.beta .* gradqy))' * gradqy;

            El.B1_beta2_loc{ie_sh} = El.B1_beta2_loc{ie_sh} + (dx .* (par.m .* par.beta.^2 .* gradqx))' * gradqx;
            El.B2_beta2_loc{ie_sh} = El.B2_beta2_loc{ie_sh} + (dx .* (par.m .* par.beta.^2 .* gradqx))' * gradqy;
            El.B3_beta2_loc{ie_sh} = El.B3_beta2_loc{ie_sh} + (dx .* (par.m .* par.beta.^2 .* gradqy))' * gradqx;
            El.B4_beta2_loc{ie_sh} = El.B4_beta2_loc{ie_sh} + (dx .* (par.m .* par.beta.^2 .* gradqy))' * gradqy;

            El.P1_loc{ie_sh} = El.P1_loc{ie_sh} - (dx .* (par.m .* phiq))' * gradqx;
            El.P2_loc{ie_sh} = El.P2_loc{ie_sh} - (dx .* (par.m .* phiq))' * gradqy;

            El.P1_beta_loc{ie_sh} = El.P1_beta_loc{ie_sh} - (dx .* (par.m .* par.beta .* phiq))' * gradqx;
            El.P2_beta_loc{ie_sh} = El.P2_beta_loc{ie_sh} - (dx .* (par.m .* par.beta .* phiq))' * gradqy;
            
        end
 
        %% Boundary integrals and stabilization terms
        
        % Computation of all the penalty coefficients for the element ie
        [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

        for iedg = 1 : neighbor.nedges(ie) % loop over faces
            
            % Extraction of tag and id of neighbor el
            id_el_neigh = neigh_ie(iedg);
            idneigh = (neigh_ie_unq == neighbor.neigh{ie}(iedg));

            % Extraction of the indexes for assembling face matrices (contribution of the element ie)
            ii_index_neigh{ie_sh}(1:femregion.nbases,:) = repmat(index, 1 ,femregion.nbases);
            jj_index_neigh{ie_sh}(1:femregion.nbases,:) = repmat(index',femregion.nbases,1);

            if id_el_neigh > 0
                id_neigh  = femregion.id(id_el_neigh);
                tag_neigh = femregion.tag(id_el_neigh);
            else
                id_neigh  = 0;
                tag_neigh = 'NaN';
            end
            
            % Extraction of the edge coordinates
            if iedg == neighbor.nedges(ie)
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(1,:);
            else
                p1 = coords_ie(iedg,:);
                p2 = coords_ie(iedg+1,:);
            end
            
            % Construction of quadrature nodes on the face
            [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);
            
            xq = qNodes_1D(:,1);
            yq = qNodes_1D(:,2);
            
            % Scaled weights
            ds = meshsize(iedg) * w_1D;
            
            % Extraction of normals to the face
            nx = normals(1,iedg);
            ny = normals(2,iedg);
           
            % Evaluation of physical parameters
            par.mu        = Data.mu{id_ie}(xq,yq);
            par.lam       = Data.lam{id_ie}(xq,yq);
            par.beta      = Data.beta{id_ie}(xq,yq);
            par.m         = Data.m{id_ie}(xq,yq);

            % Auxiliary quantities (cf. physical parameters) for vector assembling
            if strcmp(tag_neigh,'E') %elastic neighbor
                par.mu_n   = Data.mu_el{id_neigh-id_elas}(xq,yq);
                par.lam_n  = Data.lam_el{id_neigh-id_elas}(xq,yq);
                par.m_n    = Data.m{id_ie}(xq,yq);
                par.beta_n = Data.beta{id_ie}(xq,yq);
            elseif strcmp(tag_neigh,'P') %poroelastic neighbor
                par.mu_n   = Data.mu{id_neigh}(xq,yq);
                par.lam_n  = Data.lam{id_neigh}(xq,yq);
                par.m_n    = Data.m{id_neigh}(xq,yq);
                par.beta_n = Data.beta{id_neigh}(xq,yq);
            elseif strcmp(tag_neigh,'A') %acoustic edge
                par.mu_n      = Data.mu_el{id_ie}(xq,yq);
                par.lam_n     = Data.lam_el{id_ie}(xq,yq);
                par.m_n       = Data.m{id_ie}(xq,yq);
                par.beta_n    = Data.beta{id_ie}(xq,yq);
            elseif strcmp(tag_neigh,'NaN') %boundary edge
                par.mu_n      = Data.mu{id_ie}(xq,yq);
                par.lam_n     = Data.lam{id_ie}(xq,yq);
                par.m_n       = Data.m{id_ie}(xq,yq);
                par.beta_n    = Data.beta{id_ie}(xq,yq);
                par.vp_poroI  = Data.vp_poroI{id_ie}(xq,yq);
                par.vp_poroII = Data.vp_poroII{id_ie}(xq,yq);
                par.vs_poro   = Data.vs_poro{id_ie}(xq,yq);
                par.phi_por   = Data.phi_por{id_ie}(xq,yq);
                par.a_coef    = Data.a_coef{id_ie}(xq,yq);
            end
            par.lambda_ave = 2*par.lam .* par.lam_n ./ (par.lam + par.lam_n);
            par.mu_ave     = 2*par.mu .* par.mu_n ./ (par.mu + par.mu_n);
            par.harm_ave   = (par.lambda_ave + 2*par.mu_ave);
            par.m_ave      = 2*par.m .* par.m_n ./ (par.m + par.m_n);
            par.beta_ave   = 2*par.beta .* par.beta_n ./ (par.beta + par.beta_n);
                        
            % Construction and evalutation on the quadrature points of the basis functions
            [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);            
            
            % Dirichlet boundary faces
            if neigh_ie(iedg) == -1

                Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
                Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

                Edge.S1_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.S1_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.S1_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) + (ds .* ((par.lambda_ave + 2*par.mu_ave) .* nx .* gradedgeqx + par.mu_ave .* ny .* gradedgeqy ))' * phiedgeq;
                Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.lambda_ave .* ny .* gradedgeqx + par.mu_ave .* nx .* gradedgeqy ))' * phiedgeq;
                Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.mu_ave .* ny .* gradedgeqx + par.lambda_ave .* nx .* gradedgeqy ))' * phiedgeq;
                Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.mu_ave .* nx .* gradedgeqx + (par.lambda_ave + 2*par.mu_ave) .* ny .* gradedgeqy ))' * phiedgeq;

                Edge.BT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .*nx .* gradedgeqx ))' * phiedgeq;
                Edge.BT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .*ny .* gradedgeqx ))' * phiedgeq;
                Edge.BT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .*nx .* gradedgeqy ))' * phiedgeq;
                Edge.BT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .*ny .* gradedgeqy ))' * phiedgeq;

                Edge.BT1_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_beta_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave .*nx .* gradedgeqx ))' * phiedgeq;
                Edge.BT2_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_beta_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave .*ny .* gradedgeqx ))' * phiedgeq;
                Edge.BT3_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_beta_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave .*nx .* gradedgeqy ))' * phiedgeq;
                Edge.BT4_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_beta_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave .*ny .* gradedgeqy ))' * phiedgeq;

                Edge.BT1_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_beta2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave.^2 .*nx .* gradedgeqx ))' * phiedgeq;
                Edge.BT2_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_beta2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave.^2 .*ny .* gradedgeqx ))' * phiedgeq;
                Edge.BT3_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_beta2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave.^2 .*nx .* gradedgeqy ))' * phiedgeq;
                Edge.BT4_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_beta2_loc{ie_sh}(1:femregion.nbases,:) + (ds .* (par.m_ave .* par.beta_ave.^2 .*ny .* gradedgeqy ))' * phiedgeq;

            % Absorbing boundary faces
            elseif neigh_ie(iedg) == -3 

                Edge.ABC_uu_1_loc{ie_sh} = Edge.ABC_uu_1_loc{ie_sh} + (ds .* (par.rho.*par.vp_poroI*nx*nx + (par.rho - par.rho_f.*par.phi_por./par.a_coef).*par.vs_poro*(1-nx*nx)) .* phiedgeq)' * phiedgeq;
                Edge.ABC_uu_2_loc{ie_sh} = Edge.ABC_uu_2_loc{ie_sh} + (ds .* (par.rho.*par.vp_poroI*nx*ny - (par.rho - par.rho_f.*par.phi_por./par.a_coef).*par.vs_poro*nx*ny)     .* phiedgeq)' * phiedgeq;
                Edge.ABC_uu_3_loc{ie_sh} = Edge.ABC_uu_3_loc{ie_sh} + (ds .* (par.rho.*par.vp_poroI*ny*nx - (par.rho - par.rho_f.*par.phi_por./par.a_coef).*par.vs_poro*ny*nx)     .* phiedgeq)' * phiedgeq;
                Edge.ABC_uu_4_loc{ie_sh} = Edge.ABC_uu_4_loc{ie_sh} + (ds .* (par.rho.*par.vp_poroI*ny*ny + (par.rho - par.rho_f.*par.phi_por./par.a_coef).*par.vs_poro*(1-ny*ny)) .* phiedgeq)' * phiedgeq;

                Edge.ABC_uw_1_loc{ie_sh} = Edge.ABC_uu_1_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroII*nx*nx) .* phiedgeq)' * phiedgeq;
                Edge.ABC_uw_2_loc{ie_sh} = Edge.ABC_uu_2_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroII*nx*ny) .* phiedgeq)' * phiedgeq;
                Edge.ABC_uw_3_loc{ie_sh} = Edge.ABC_uu_3_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroII*ny*nx) .* phiedgeq)' * phiedgeq;
                Edge.ABC_uw_4_loc{ie_sh} = Edge.ABC_uu_4_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroII*ny*ny) .* phiedgeq)' * phiedgeq;

                Edge.ABC_wu_1_loc{ie_sh} = Edge.ABC_uu_1_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroI*nx*nx) .* phiedgeq)' * phiedgeq;
                Edge.ABC_wu_2_loc{ie_sh} = Edge.ABC_uu_2_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroI*nx*ny) .* phiedgeq)' * phiedgeq;
                Edge.ABC_wu_3_loc{ie_sh} = Edge.ABC_uu_3_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroI*ny*nx) .* phiedgeq)' * phiedgeq;
                Edge.ABC_wu_4_loc{ie_sh} = Edge.ABC_uu_4_loc{ie_sh} + (ds .* (par.rho_f.*par.vp_poroI*ny*ny) .* phiedgeq)' * phiedgeq;

                Edge.ABC_ww_1_loc{ie_sh} = Edge.ABC_uu_1_loc{ie_sh} + (ds .* (par.rho_w.*par.vp_poroII*nx*nx) .* phiedgeq)' * phiedgeq;
                Edge.ABC_ww_2_loc{ie_sh} = Edge.ABC_uu_2_loc{ie_sh} + (ds .* (par.rho_w.*par.vp_poroII*nx*ny) .* phiedgeq)' * phiedgeq;
                Edge.ABC_ww_3_loc{ie_sh} = Edge.ABC_uu_3_loc{ie_sh} + (ds .* (par.rho_w.*par.vp_poroII*ny*nx) .* phiedgeq)' * phiedgeq;
                Edge.ABC_ww_4_loc{ie_sh} = Edge.ABC_uu_4_loc{ie_sh} + (ds .* (par.rho_w.*par.vp_poroII*ny*ny) .* phiedgeq)' * phiedgeq;
                                             
            % Poroelastic neighbor
            elseif neigh_ie(iedg)>0 && femregion.tag(neigh_ie(iedg)) == 'P'
                
                % Element itself
                Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
                Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_P_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

                Edge.S1_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .*nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .*nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .*ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .*ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.S1_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.S1_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT1_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* ((par.lambda_ave + 2*par.mu_ave) .* nx .* gradedgeqx + par.mu_ave .* ny .* gradedgeqy ))' * phiedgeq;
                Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.lambda_ave .* ny .* gradedgeqx + par.mu_ave .* nx .* gradedgeqy ))' * phiedgeq;
                Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT3_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.mu_ave .* ny .* gradedgeqx + par.lambda_ave .* nx .* gradedgeqy ))' * phiedgeq;
                Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.IT4_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.mu_ave .* nx .* gradedgeqx + (par.lambda_ave + 2*par.mu_ave) .* ny .* gradedgeqy ))' * phiedgeq;

                Edge.BT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .*nx .* gradedgeqx ))' * phiedgeq;
                Edge.BT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .*ny .* gradedgeqx ))' * phiedgeq;
                Edge.BT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .*nx .* gradedgeqy ))' * phiedgeq;
                Edge.BT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .*ny .* gradedgeqy ))' * phiedgeq;

                Edge.BT1_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave .*nx .* gradedgeqx ))' * phiedgeq;
                Edge.BT2_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave .*ny .* gradedgeqx ))' * phiedgeq;
                Edge.BT3_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave .*nx .* gradedgeqy ))' * phiedgeq;
                Edge.BT4_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave .*ny .* gradedgeqy ))' * phiedgeq;

                Edge.BT1_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*nx .* gradedgeqx ))' * phiedgeq;
                Edge.BT2_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*ny .* gradedgeqx ))' * phiedgeq;
                Edge.BT3_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*nx .* gradedgeqy ))' * phiedgeq;
                Edge.BT4_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*ny .* gradedgeqy ))' * phiedgeq;

                % Neighboring element
                neigh_idx = find(idneigh)*femregion.nbases+1:(find(idneigh)+1)*(femregion.nbases);
                index_neigh = (neighbor.neigh{ie}(iedg)-1-nel_sh)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

                % Extraction of the indexes for assembling face matrices (contribution of the neighboring element)
                ii_index_neigh{ie_sh}(neigh_idx,:) = repmat(index, 1,femregion.nbases);
                jj_index_neigh{ie_sh}(neigh_idx,:) = repmat(index_neigh',femregion.nbases,1);

                % Construction and evalutation on the quadrature points of the basis functions for the neighbor
                phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);
                
                % Neighboring element
                Edge.S1_P_loc{ie_sh}(neigh_idx,:) = Edge.S1_P_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeqneigh;
                Edge.S4_P_loc{ie_sh}(neigh_idx,:) = Edge.S4_P_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeqneigh;

                Edge.S1_B_loc{ie_sh}(neigh_idx,:) = Edge.S1_B_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .*nx .* nx .* phiedgeq)' * phiedgeqneigh;
                Edge.S2_B_loc{ie_sh}(neigh_idx,:) = Edge.S2_B_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .*nx .* ny .* phiedgeq)' * phiedgeqneigh;
                Edge.S3_B_loc{ie_sh}(neigh_idx,:) = Edge.S3_B_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .*ny .* nx .* phiedgeq)' * phiedgeqneigh;
                Edge.S4_B_loc{ie_sh}(neigh_idx,:) = Edge.S4_B_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .*ny .* ny .* phiedgeq)' * phiedgeqneigh;

                Edge.S1_B_beta_loc{ie_sh}(neigh_idx,:) = Edge.S1_B_beta_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*nx .* nx .* phiedgeq)' * phiedgeqneigh;
                Edge.S2_B_beta_loc{ie_sh}(neigh_idx,:) = Edge.S2_B_beta_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*nx .* ny .* phiedgeq)' * phiedgeqneigh;
                Edge.S3_B_beta_loc{ie_sh}(neigh_idx,:) = Edge.S3_B_beta_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*ny .* nx .* phiedgeq)' * phiedgeqneigh;
                Edge.S4_B_beta_loc{ie_sh}(neigh_idx,:) = Edge.S4_B_beta_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*ny .* ny .* phiedgeq)' * phiedgeqneigh;

                Edge.S1_B_beta2_loc{ie_sh}(neigh_idx,:) = Edge.S1_B_beta2_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*nx .* nx .* phiedgeq)' * phiedgeqneigh;
                Edge.S2_B_beta2_loc{ie_sh}(neigh_idx,:) = Edge.S2_B_beta2_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*nx .* ny .* phiedgeq)' * phiedgeqneigh;
                Edge.S3_B_beta2_loc{ie_sh}(neigh_idx,:) = Edge.S3_B_beta2_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*ny .* nx .* phiedgeq)' * phiedgeqneigh;
                Edge.S4_B_beta2_loc{ie_sh}(neigh_idx,:) = Edge.S4_B_beta2_loc{ie_sh}(neigh_idx,:) - penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*ny .* ny .* phiedgeq)' * phiedgeqneigh;

                Edge.IT1_loc{ie_sh}(neigh_idx,:) = Edge.IT1_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* ((par.lambda_ave + 2*par.mu_ave) .* nx .* gradedgeqx + par.mu_ave .* ny .* gradedgeqy ))' * phiedgeqneigh;
                Edge.IT2_loc{ie_sh}(neigh_idx,:) = Edge.IT2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.lambda_ave .* ny .* gradedgeqx + par.mu_ave .* nx .* gradedgeqy ))' * phiedgeqneigh;
                Edge.IT3_loc{ie_sh}(neigh_idx,:) = Edge.IT3_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.mu_ave .* ny .* gradedgeqx + par.lambda_ave .* nx .* gradedgeqy ))' * phiedgeqneigh;
                Edge.IT4_loc{ie_sh}(neigh_idx,:) = Edge.IT4_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.mu_ave .* nx .* gradedgeqx + (par.lambda_ave + 2*par.mu_ave) .* ny .* gradedgeqy ))' * phiedgeqneigh;

                Edge.BT1_loc{ie_sh}(neigh_idx,:) = Edge.BT1_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .*nx .* gradedgeqx ))' * phiedgeqneigh;
                Edge.BT2_loc{ie_sh}(neigh_idx,:) = Edge.BT2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .*ny .* gradedgeqx ))' * phiedgeqneigh;
                Edge.BT3_loc{ie_sh}(neigh_idx,:) = Edge.BT3_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .*nx .* gradedgeqy ))' * phiedgeqneigh;
                Edge.BT4_loc{ie_sh}(neigh_idx,:) = Edge.BT4_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .*ny .* gradedgeqy ))' * phiedgeqneigh;

                Edge.BT1_beta_loc{ie_sh}(neigh_idx,:) = Edge.BT1_beta_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave .*nx .* gradedgeqx ))' * phiedgeqneigh;
                Edge.BT2_beta_loc{ie_sh}(neigh_idx,:) = Edge.BT2_beta_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave .*ny .* gradedgeqx ))' * phiedgeqneigh;
                Edge.BT3_beta_loc{ie_sh}(neigh_idx,:) = Edge.BT3_beta_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave .*nx .* gradedgeqy ))' * phiedgeqneigh;
                Edge.BT4_beta_loc{ie_sh}(neigh_idx,:) = Edge.BT4_beta_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave .*ny .* gradedgeqy ))' * phiedgeqneigh;

                Edge.BT1_beta2_loc{ie_sh}(neigh_idx,:) = Edge.BT1_beta2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*nx .* gradedgeqx ))' * phiedgeqneigh;
                Edge.BT2_beta2_loc{ie_sh}(neigh_idx,:) = Edge.BT2_beta2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*ny .* gradedgeqx ))' * phiedgeqneigh;
                Edge.BT3_beta2_loc{ie_sh}(neigh_idx,:) = Edge.BT3_beta2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*nx .* gradedgeqy ))' * phiedgeqneigh;
                Edge.BT4_beta2_loc{ie_sh}(neigh_idx,:) = Edge.BT4_beta2_loc{ie_sh}(neigh_idx,:) - 0.5 * (ds .* (par.m_ave .* par.beta_ave.^2 .*ny .* gradedgeqy ))' * phiedgeqneigh;
                
            % Elastic neighbor
            elseif neigh_ie(iedg)>0 && femregion.tag(neigh_ie(iedg)) == 'E'
                
                % Element itself
                Edge.S1_P_EP_loc{ie_sh} = Edge.S1_P_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
                Edge.S4_P_EP_loc{ie_sh} = Edge.S4_P_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

                Edge.BT1_beta2_EP_loc{ie_sh} = Edge.BT1_beta2_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* phiedgeq)' * gradedgeqx;
                Edge.BT2_beta2_EP_loc{ie_sh} = Edge.BT2_beta2_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* phiedgeq)' * gradedgeqy;
                Edge.BT3_beta2_EP_loc{ie_sh} = Edge.BT3_beta2_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* phiedgeq)' * gradedgeqx;
                Edge.BT4_beta2_EP_loc{ie_sh} = Edge.BT4_beta2_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* phiedgeq)' * gradedgeqy;

                Edge.S1_B_beta2_EP_loc{ie_sh} = Edge.S1_B_beta2_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_beta2_EP_loc{ie_sh} = Edge.S2_B_beta2_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_beta2_EP_loc{ie_sh} = Edge.S3_B_beta2_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_beta2_EP_loc{ie_sh} = Edge.S4_B_beta2_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .*ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.S1_B_beta_EP_loc{ie_sh} = Edge.S1_B_beta_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_beta_EP_loc{ie_sh} = Edge.S2_B_beta_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_beta_EP_loc{ie_sh} = Edge.S3_B_beta_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_beta_EP_loc{ie_sh} = Edge.S4_B_beta_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .*ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.BT1_beta_EP_loc{ie_sh} = Edge.BT1_beta_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave .* nx .* gradedgeqx)' *  phiedgeq;
                Edge.BT2_beta_EP_loc{ie_sh} = Edge.BT2_beta_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave .* ny .* gradedgeqx)' *  phiedgeq;
                Edge.BT3_beta_EP_loc{ie_sh} = Edge.BT3_beta_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave .* nx .* gradedgeqy)' *  phiedgeq;
                Edge.BT4_beta_EP_loc{ie_sh} = Edge.BT4_beta_EP_loc{ie_sh} + (ds .* par.m_ave .* par.beta_ave .* ny .* gradedgeqy)' *  phiedgeq;

                Edge.S1_B_EP_loc{ie_sh} = Edge.S1_B_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .*nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.S2_B_EP_loc{ie_sh} = Edge.S2_B_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .*nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.S3_B_EP_loc{ie_sh} = Edge.S3_B_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .*ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.S4_B_EP_loc{ie_sh} = Edge.S4_B_EP_loc{ie_sh} + penalty_geom(iedg) * (ds .* par.m_ave .*ny .* ny .* phiedgeq)' * phiedgeq;

                Edge.BT1_EP_loc{ie_sh} = Edge.BT1_EP_loc{ie_sh} + (ds .* par.m_ave .* nx .* gradedgeqx)' *  phiedgeq;
                Edge.BT2_EP_loc{ie_sh} = Edge.BT2_EP_loc{ie_sh} + (ds .* par.m_ave .* ny .* gradedgeqx)' *  phiedgeq;
                Edge.BT3_EP_loc{ie_sh} = Edge.BT3_EP_loc{ie_sh} + (ds .* par.m_ave .* nx .* gradedgeqy)' *  phiedgeq;
                Edge.BT4_EP_loc{ie_sh} = Edge.BT4_EP_loc{ie_sh} + (ds .* par.m_ave .* ny .* gradedgeqy)' *  phiedgeq;

            % Acoustic neighbor
            elseif neigh_ie(iedg) >0 && femregion.tag(neigh_ie(iedg)) == 'A'

                % Element itself
                Edge.D1_loc{ie_sh} = Edge.D1_loc{ie_sh} + (ds .* nx .* nx .* phiedgeq)' * phiedgeq;
                Edge.D2_loc{ie_sh} = Edge.D2_loc{ie_sh} + (ds .* nx .* ny .* phiedgeq)' * phiedgeq;
                Edge.D3_loc{ie_sh} = Edge.D3_loc{ie_sh} + (ds .* ny .* nx .* phiedgeq)' * phiedgeq;
                Edge.D4_loc{ie_sh} = Edge.D4_loc{ie_sh} + (ds .* ny .* ny .* phiedgeq)' * phiedgeq;


                if (Data.tau == 0)

                    Edge.S1_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* nx .* nx .* phiedgeq)' * phiedgeq;
                    Edge.S2_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* nx .* ny .* phiedgeq)' * phiedgeq;
                    Edge.S3_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* ny .* nx .* phiedgeq)' * phiedgeq;
                    Edge.S4_B_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* ny .* ny .* phiedgeq)' * phiedgeq;

                    Edge.S1_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* nx .* nx .* phiedgeq)' * phiedgeq;
                    Edge.S2_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* nx .* ny .* phiedgeq)' * phiedgeq;
                    Edge.S3_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* ny .* nx .* phiedgeq)' * phiedgeq;
                    Edge.S4_B_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_beta_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave .* ny .* ny .* phiedgeq)' * phiedgeq;

                    Edge.S1_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S1_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* nx .* phiedgeq)' * phiedgeq;
                    Edge.S2_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S2_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* ny .* phiedgeq)' * phiedgeq;
                    Edge.S3_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S3_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* nx .* phiedgeq)' * phiedgeq;
                    Edge.S4_B_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.S4_B_beta2_loc{ie_sh}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* ny .* phiedgeq)' * phiedgeq;

                    Edge.BT1_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* nx .* gradedgeqx)' * phiedgeq;
                    Edge.BT2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* ny .* gradedgeqx)' * phiedgeq;
                    Edge.BT3_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* nx .* gradedgeqy)' * phiedgeq;
                    Edge.BT4_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* ny .* gradedgeqy)' * phiedgeq;

                    Edge.BT1_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave .* nx .* gradedgeqx)' * phiedgeq;
                    Edge.BT2_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave .* ny .* gradedgeqx)' * phiedgeq;
                    Edge.BT3_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave .* nx .* gradedgeqy)' * phiedgeq;
                    Edge.BT4_beta_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_beta_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave .* ny .* gradedgeqy)' * phiedgeq;

                    Edge.BT1_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT1_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* gradedgeqx)' * phiedgeq;
                    Edge.BT2_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT2_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* gradedgeqx)' * phiedgeq;
                    Edge.BT3_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT3_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave.^2 .* nx .* gradedgeqy)' * phiedgeq;
                    Edge.BT4_beta2_loc{ie_sh}(1:femregion.nbases,:) = Edge.BT4_beta2_loc{ie_sh}(1:femregion.nbases,:) + 0.5 * (ds .* par.m_ave .* par.beta_ave.^2 .* ny .* gradedgeqy)' * phiedgeq;

                end
            end

        end

    end
    
end

% Local matrix to global matrix

%% Reshape of the volume matrices 
ii_index  = reshape(cell2mat(ii_index),[femregion.nbases,femregion.nbases*femregion.nel_p]);
jj_index  = reshape(cell2mat(jj_index),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.V1_loc = reshape(cell2mat(El.V1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.V2_loc = reshape(cell2mat(El.V2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.V3_loc = reshape(cell2mat(El.V3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.V4_loc = reshape(cell2mat(El.V4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.M1_P_rhof_loc = reshape(cell2mat(El.M1_P_rhof_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.M1_P_rhow_loc = reshape(cell2mat(El.M1_P_rhow_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.M1_P_rho_loc  = reshape(cell2mat(El.M1_P_rho_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.M1_P_eta_kper_loc  = reshape(cell2mat(El.M1_P_eta_kper_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.M1_P_rho2_zeta_loc = reshape(cell2mat(El.M1_P_rho2_zeta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.M1_P_rho_zeta2_loc = reshape(cell2mat(El.M1_P_rho_zeta2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.MPrjP_1_loc   = reshape(cell2mat(El.MPrjP_1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.B1_loc = reshape(cell2mat(El.B1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B2_loc = reshape(cell2mat(El.B2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B3_loc = reshape(cell2mat(El.B3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B4_loc = reshape(cell2mat(El.B4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.B1_beta_loc = reshape(cell2mat(El.B1_beta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B2_beta_loc = reshape(cell2mat(El.B2_beta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B3_beta_loc = reshape(cell2mat(El.B3_beta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B4_beta_loc = reshape(cell2mat(El.B4_beta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.B1_beta2_loc = reshape(cell2mat(El.B1_beta2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B2_beta2_loc = reshape(cell2mat(El.B2_beta2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B3_beta2_loc = reshape(cell2mat(El.B3_beta2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.B4_beta2_loc = reshape(cell2mat(El.B4_beta2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.P1_loc = reshape(cell2mat(El.P1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.P2_loc = reshape(cell2mat(El.P2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

El.P1_beta_loc = reshape(cell2mat(El.P1_beta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
El.P2_beta_loc = reshape(cell2mat(El.P2_beta_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

%% Reshape of the edge matrices 

ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.S1_P_loc  = reshape(cell2mat(Edge.S1_P_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S4_P_loc  = reshape(cell2mat(Edge.S4_P_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.IT1_loc  = reshape(cell2mat(Edge.IT1_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.IT2_loc  = reshape(cell2mat(Edge.IT2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.IT3_loc  = reshape(cell2mat(Edge.IT3_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.IT4_loc  = reshape(cell2mat(Edge.IT4_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.S1_B_loc = reshape(cell2mat(Edge.S1_B_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S2_B_loc = reshape(cell2mat(Edge.S2_B_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S3_B_loc = reshape(cell2mat(Edge.S3_B_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S4_B_loc = reshape(cell2mat(Edge.S4_B_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.S1_B_beta_loc = reshape(cell2mat(Edge.S1_B_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S2_B_beta_loc = reshape(cell2mat(Edge.S2_B_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S3_B_beta_loc = reshape(cell2mat(Edge.S3_B_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S4_B_beta_loc = reshape(cell2mat(Edge.S4_B_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.S1_B_beta2_loc = reshape(cell2mat(Edge.S1_B_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S2_B_beta2_loc = reshape(cell2mat(Edge.S2_B_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S3_B_beta2_loc = reshape(cell2mat(Edge.S3_B_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.S4_B_beta2_loc = reshape(cell2mat(Edge.S4_B_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.BT1_loc  = reshape(cell2mat(Edge.BT1_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT2_loc  = reshape(cell2mat(Edge.BT2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT3_loc  = reshape(cell2mat(Edge.BT3_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT4_loc  = reshape(cell2mat(Edge.BT4_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.BT1_beta_loc  = reshape(cell2mat(Edge.BT1_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT2_beta_loc  = reshape(cell2mat(Edge.BT2_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT3_beta_loc  = reshape(cell2mat(Edge.BT3_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT4_beta_loc  = reshape(cell2mat(Edge.BT4_beta_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

Edge.BT1_beta2_loc  = reshape(cell2mat(Edge.BT1_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT2_beta2_loc  = reshape(cell2mat(Edge.BT2_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT3_beta2_loc  = reshape(cell2mat(Edge.BT3_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);
Edge.BT4_beta2_loc  = reshape(cell2mat(Edge.BT4_beta2_loc),[femregion.nbases,femregion.nel_p*max_nedges*femregion.nbases]);

%% Reshape of the absorbing boundary matrices

Edge.ABC_uu_1_loc = reshape(cell2mat(Edge.ABC_uu_1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_uu_2_loc = reshape(cell2mat(Edge.ABC_uu_2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_uu_3_loc = reshape(cell2mat(Edge.ABC_uu_3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_uu_4_loc = reshape(cell2mat(Edge.ABC_uu_4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.ABC_wu_1_loc = reshape(cell2mat(Edge.ABC_wu_1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_wu_2_loc = reshape(cell2mat(Edge.ABC_wu_2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_wu_3_loc = reshape(cell2mat(Edge.ABC_wu_3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_wu_4_loc = reshape(cell2mat(Edge.ABC_wu_4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.ABC_uw_1_loc = reshape(cell2mat(Edge.ABC_uw_1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_uw_2_loc = reshape(cell2mat(Edge.ABC_uw_2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_uw_3_loc = reshape(cell2mat(Edge.ABC_uw_3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_uw_4_loc = reshape(cell2mat(Edge.ABC_uw_4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.ABC_ww_1_loc = reshape(cell2mat(Edge.ABC_ww_1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_ww_2_loc = reshape(cell2mat(Edge.ABC_ww_2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_ww_3_loc = reshape(cell2mat(Edge.ABC_ww_3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.ABC_ww_4_loc = reshape(cell2mat(Edge.ABC_ww_4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

%% Reshape of the acoustic boundary matrices

Edge.D1_loc = reshape(cell2mat(Edge.D1_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.D2_loc = reshape(cell2mat(Edge.D2_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.D3_loc = reshape(cell2mat(Edge.D3_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.D4_loc = reshape(cell2mat(Edge.D4_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
 
%% Reshape of the elastic boundary matrices

Edge.S1_P_EP_loc = reshape(cell2mat(Edge.S1_P_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S4_P_EP_loc = reshape(cell2mat(Edge.S4_P_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
 
Edge.S1_B_EP_loc = reshape(cell2mat(Edge.S1_B_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S2_B_EP_loc = reshape(cell2mat(Edge.S2_B_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S3_B_EP_loc = reshape(cell2mat(Edge.S3_B_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S4_B_EP_loc = reshape(cell2mat(Edge.S4_B_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.S1_B_beta_EP_loc = reshape(cell2mat(Edge.S1_B_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S2_B_beta_EP_loc = reshape(cell2mat(Edge.S2_B_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S3_B_beta_EP_loc = reshape(cell2mat(Edge.S3_B_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S4_B_beta_EP_loc = reshape(cell2mat(Edge.S4_B_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.S1_B_beta2_EP_loc = reshape(cell2mat(Edge.S1_B_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S2_B_beta2_EP_loc = reshape(cell2mat(Edge.S2_B_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S3_B_beta2_EP_loc = reshape(cell2mat(Edge.S3_B_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.S4_B_beta2_EP_loc = reshape(cell2mat(Edge.S4_B_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.BT1_EP_loc = reshape(cell2mat(Edge.BT1_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT2_EP_loc = reshape(cell2mat(Edge.BT2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT3_EP_loc = reshape(cell2mat(Edge.BT3_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT4_EP_loc = reshape(cell2mat(Edge.BT4_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.BT1_beta_EP_loc = reshape(cell2mat(Edge.BT1_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT2_beta_EP_loc = reshape(cell2mat(Edge.BT2_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT3_beta_EP_loc = reshape(cell2mat(Edge.BT3_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT4_beta_EP_loc = reshape(cell2mat(Edge.BT4_beta_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);

Edge.BT1_beta2_EP_loc = reshape(cell2mat(Edge.BT1_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT2_beta2_EP_loc = reshape(cell2mat(Edge.BT2_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT3_beta2_EP_loc = reshape(cell2mat(Edge.BT3_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);
Edge.BT4_beta2_EP_loc = reshape(cell2mat(Edge.BT4_beta2_EP_loc),[femregion.nbases,femregion.nbases*femregion.nel_p]);


%% Elimination of additional rows in edge matrices

del = all(jj_index_neigh == 0,1);
ii_index_neigh(:,del) = [];
jj_index_neigh(:,del) = [];

Edge.S1_P_loc(:,del) = [];
Edge.S4_P_loc(:,del) = [];

Edge.IT1_loc(:,del) = [];
Edge.IT2_loc(:,del) = [];
Edge.IT3_loc(:,del) = [];
Edge.IT4_loc(:,del) = [];

Edge.S1_B_loc(:,del) = [];
Edge.S2_B_loc(:,del) = [];
Edge.S3_B_loc(:,del) = [];
Edge.S4_B_loc(:,del) = [];

Edge.S1_B_beta_loc(:,del) = [];
Edge.S2_B_beta_loc(:,del) = [];
Edge.S3_B_beta_loc(:,del) = [];
Edge.S4_B_beta_loc(:,del) = [];

Edge.S1_B_beta2_loc(:,del) = [];
Edge.S2_B_beta2_loc(:,del) = [];
Edge.S3_B_beta2_loc(:,del) = [];
Edge.S4_B_beta2_loc(:,del) = [];

Edge.BT1_loc(:,del) = [];
Edge.BT2_loc(:,del) = [];
Edge.BT3_loc(:,del) = [];
Edge.BT4_loc(:,del) = [];

Edge.BT1_beta_loc(:,del) = [];
Edge.BT2_beta_loc(:,del) = [];
Edge.BT3_beta_loc(:,del) = [];
Edge.BT4_beta_loc(:,del) = [];

Edge.BT1_beta2_loc(:,del) = [];
Edge.BT2_beta2_loc(:,del) = [];
Edge.BT3_beta2_loc(:,del) = [];
Edge.BT4_beta2_loc(:,del) = [];

%% Creation of global volume matrices 

A.V1 = sparse(ii_index,jj_index,El.V1_loc,femregion.ndof_p,femregion.ndof_p);
A.V2 = sparse(ii_index,jj_index,El.V2_loc,femregion.ndof_p,femregion.ndof_p);
A.V3 = sparse(ii_index,jj_index,El.V3_loc,femregion.ndof_p,femregion.ndof_p);
A.V4 = sparse(ii_index,jj_index,El.V4_loc,femregion.ndof_p,femregion.ndof_p);

A.M1_P_rhof  = sparse(ii_index,jj_index,El.M1_P_rhof_loc,femregion.ndof_p,femregion.ndof_p);
A.M1_P_rhow  = sparse(ii_index,jj_index,El.M1_P_rhow_loc,femregion.ndof_p,femregion.ndof_p);
A.M1_P_rho  = sparse(ii_index,jj_index,El.M1_P_rho_loc,femregion.ndof_p,femregion.ndof_p);

A.M1_P_eta_kper  = sparse(ii_index,jj_index,El.M1_P_eta_kper_loc,femregion.ndof_p,femregion.ndof_p);
A.M1_P_rho2_zeta = sparse(ii_index,jj_index,El.M1_P_rho2_zeta_loc,femregion.ndof_p,femregion.ndof_p);
A.M1_P_rho_zeta2 = sparse(ii_index,jj_index,El.M1_P_rho_zeta2_loc,femregion.ndof_p,femregion.ndof_p);

A.MPrjP_1  = sparse(ii_index,jj_index,El.MPrjP_1_loc,femregion.ndof_p,femregion.ndof_p);

A.B1 = sparse(ii_index,jj_index,El.B1_loc,femregion.ndof_p,femregion.ndof_p);
A.B2 = sparse(ii_index,jj_index,El.B2_loc,femregion.ndof_p,femregion.ndof_p);
A.B3 = sparse(ii_index,jj_index,El.B3_loc,femregion.ndof_p,femregion.ndof_p);
A.B4 = sparse(ii_index,jj_index,El.B4_loc,femregion.ndof_p,femregion.ndof_p);

A.B1_beta = sparse(ii_index,jj_index,El.B1_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.B2_beta = sparse(ii_index,jj_index,El.B2_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.B3_beta = sparse(ii_index,jj_index,El.B3_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.B4_beta = sparse(ii_index,jj_index,El.B4_beta_loc,femregion.ndof_p,femregion.ndof_p);

A.B1_beta2 = sparse(ii_index,jj_index,El.B1_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.B2_beta2 = sparse(ii_index,jj_index,El.B2_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.B3_beta2 = sparse(ii_index,jj_index,El.B3_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.B4_beta2 = sparse(ii_index,jj_index,El.B4_beta2_loc,femregion.ndof_p,femregion.ndof_p);

A.P1 = sparse(ii_index,jj_index,El.P1_loc,femregion.ndof_p,femregion.ndof_p);
A.P2 = sparse(ii_index,jj_index,El.P2_loc,femregion.ndof_p,femregion.ndof_p);

A.P1_beta = sparse(ii_index,jj_index,El.P1_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.P2_beta = sparse(ii_index,jj_index,El.P2_beta_loc,femregion.ndof_p,femregion.ndof_p);
               
%% Creation of global edge matrices 

A.S1_P     = sparse(ii_index_neigh,jj_index_neigh,Edge.S1_P_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_P     = sparse(femregion.ndof_p,femregion.ndof_p);
A.S3_P     = sparse(femregion.ndof_p,femregion.ndof_p);
A.S4_P     = sparse(ii_index_neigh,jj_index_neigh,Edge.S4_P_loc,femregion.ndof_p,femregion.ndof_p);

A.IT1_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT1_loc,femregion.ndof_p,femregion.ndof_p);
A.IT2_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT2_loc,femregion.ndof_p,femregion.ndof_p);
A.IT3_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT3_loc,femregion.ndof_p,femregion.ndof_p);
A.IT4_P    = sparse(ii_index_neigh,jj_index_neigh,Edge.IT4_loc,femregion.ndof_p,femregion.ndof_p);

A.S1_B     = sparse(ii_index_neigh,jj_index_neigh,Edge.S1_B_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_B     = sparse(ii_index_neigh,jj_index_neigh,Edge.S2_B_loc,femregion.ndof_p,femregion.ndof_p);
A.S3_B     = sparse(ii_index_neigh,jj_index_neigh,Edge.S3_B_loc,femregion.ndof_p,femregion.ndof_p);
A.S4_B     = sparse(ii_index_neigh,jj_index_neigh,Edge.S4_B_loc,femregion.ndof_p,femregion.ndof_p);

A.S1_B_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.S1_B_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_B_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.S2_B_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.S3_B_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.S3_B_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.S4_B_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.S4_B_beta_loc,femregion.ndof_p,femregion.ndof_p);

A.S1_B_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.S1_B_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_B_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.S2_B_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.S3_B_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.S3_B_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.S4_B_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.S4_B_beta2_loc,femregion.ndof_p,femregion.ndof_p);

A.BT1    = sparse(ii_index_neigh,jj_index_neigh,Edge.BT1_loc,femregion.ndof_p,femregion.ndof_p);
A.BT2    = sparse(ii_index_neigh,jj_index_neigh,Edge.BT2_loc,femregion.ndof_p,femregion.ndof_p);
A.BT3    = sparse(ii_index_neigh,jj_index_neigh,Edge.BT3_loc,femregion.ndof_p,femregion.ndof_p);
A.BT4    = sparse(ii_index_neigh,jj_index_neigh,Edge.BT4_loc,femregion.ndof_p,femregion.ndof_p);

A.BT1_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.BT1_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.BT2_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.BT2_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.BT3_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.BT3_beta_loc,femregion.ndof_p,femregion.ndof_p);
A.BT4_beta = sparse(ii_index_neigh,jj_index_neigh,Edge.BT4_beta_loc,femregion.ndof_p,femregion.ndof_p);

A.BT1_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.BT1_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.BT2_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.BT2_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.BT3_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.BT3_beta2_loc,femregion.ndof_p,femregion.ndof_p);
A.BT4_beta2 = sparse(ii_index_neigh,jj_index_neigh,Edge.BT4_beta2_loc,femregion.ndof_p,femregion.ndof_p);

%% Reshape of the absorbing boundary matrices

A.ABC_uu_1   = sparse(ii_index,jj_index,Edge.ABC_uu_1_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_uu_2   = sparse(ii_index,jj_index,Edge.ABC_uu_2_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_uu_3   = sparse(ii_index,jj_index,Edge.ABC_uu_3_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_uu_4   = sparse(ii_index,jj_index,Edge.ABC_uu_4_loc,femregion.ndof_p,femregion.ndof_p);

A.ABC_wu_1   = sparse(ii_index,jj_index,Edge.ABC_wu_1_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_wu_2   = sparse(ii_index,jj_index,Edge.ABC_wu_2_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_wu_3   = sparse(ii_index,jj_index,Edge.ABC_wu_3_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_wu_4   = sparse(ii_index,jj_index,Edge.ABC_wu_4_loc,femregion.ndof_p,femregion.ndof_p);

A.ABC_uw_1   = sparse(ii_index,jj_index,Edge.ABC_uw_1_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_uw_2   = sparse(ii_index,jj_index,Edge.ABC_uw_2_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_uw_3   = sparse(ii_index,jj_index,Edge.ABC_uw_3_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_uw_4   = sparse(ii_index,jj_index,Edge.ABC_uw_4_loc,femregion.ndof_p,femregion.ndof_p);

A.ABC_ww_1   = sparse(ii_index,jj_index,Edge.ABC_ww_1_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_ww_2   = sparse(ii_index,jj_index,Edge.ABC_ww_2_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_ww_3   = sparse(ii_index,jj_index,Edge.ABC_ww_3_loc,femregion.ndof_p,femregion.ndof_p);
A.ABC_ww_4   = sparse(ii_index,jj_index,Edge.ABC_ww_4_loc,femregion.ndof_p,femregion.ndof_p);

%% Reshape of the acoustic boundary matrices

A.D1 = sparse(ii_index,jj_index,Edge.D1_loc,femregion.ndof_p,femregion.ndof_p);
A.D2 = sparse(ii_index,jj_index,Edge.D2_loc,femregion.ndof_p,femregion.ndof_p);
A.D3 = sparse(ii_index,jj_index,Edge.D3_loc,femregion.ndof_p,femregion.ndof_p);
A.D4 = sparse(ii_index,jj_index,Edge.D4_loc,femregion.ndof_p,femregion.ndof_p);
 
%% Reshape of the elastic boundary matrices

A.S1_P_EP = sparse(ii_index,jj_index,Edge.S1_P_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_P_EP = sparse(femregion.ndof_p,femregion.ndof_p);
A.S3_P_EP = sparse(femregion.ndof_p,femregion.ndof_p);
A.S4_P_EP = sparse(ii_index,jj_index,Edge.S4_P_EP_loc,femregion.ndof_p,femregion.ndof_p);

A.S1_B_EP = sparse(ii_index,jj_index,Edge.S1_B_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_B_EP = sparse(ii_index,jj_index,Edge.S2_B_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S3_B_EP = sparse(ii_index,jj_index,Edge.S3_B_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S4_B_EP = sparse(ii_index,jj_index,Edge.S4_B_EP_loc,femregion.ndof_p,femregion.ndof_p);
 
A.S1_B_beta_EP = sparse(ii_index,jj_index,Edge.S1_B_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_B_beta_EP = sparse(ii_index,jj_index,Edge.S2_B_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S3_B_beta_EP = sparse(ii_index,jj_index,Edge.S3_B_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S4_B_beta_EP = sparse(ii_index,jj_index,Edge.S4_B_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);

A.S1_B_beta2_EP = sparse(ii_index,jj_index,Edge.S1_B_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S2_B_beta2_EP = sparse(ii_index,jj_index,Edge.S2_B_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S3_B_beta2_EP = sparse(ii_index,jj_index,Edge.S3_B_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.S4_B_beta2_EP = sparse(ii_index,jj_index,Edge.S4_B_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);

A.BT1_EP = sparse(ii_index,jj_index,Edge.BT1_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT2_EP = sparse(ii_index,jj_index,Edge.BT2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT3_EP = sparse(ii_index,jj_index,Edge.BT3_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT4_EP = sparse(ii_index,jj_index,Edge.BT4_EP_loc,femregion.ndof_p,femregion.ndof_p);

A.BT1_beta_EP = sparse(ii_index,jj_index,Edge.BT1_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT2_beta_EP = sparse(ii_index,jj_index,Edge.BT2_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT3_beta_EP = sparse(ii_index,jj_index,Edge.BT3_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT4_beta_EP = sparse(ii_index,jj_index,Edge.BT4_beta_EP_loc,femregion.ndof_p,femregion.ndof_p);

A.BT1_beta2_EP = sparse(ii_index,jj_index,Edge.BT1_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT2_beta2_EP = sparse(ii_index,jj_index,Edge.BT2_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT3_beta2_EP = sparse(ii_index,jj_index,Edge.BT3_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);
A.BT4_beta2_EP = sparse(ii_index,jj_index,Edge.BT4_beta2_EP_loc,femregion.ndof_p,femregion.ndof_p);

%% Building the matrices

%Mass matrix poroelastic
M_P_rho   = [A.M1_P_rho, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.M1_P_rho];
%
M_P_rhof  = [A.M1_P_rhof, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.M1_P_rhof];
%
M_P_rhow  = [A.M1_P_rhow, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.M1_P_rhow];
%
M_P_eta_kper  = [A.M1_P_eta_kper, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.M1_P_eta_kper];
%
D_vel = [A.M1_P_rho2_zeta, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.M1_P_rho2_zeta];
%
D_dis = [A.M1_P_rho_zeta2, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.M1_P_rho_zeta2];

%Projection matrix-poroelastic: here equal to mass matrix
MprjP   = [A.MPrjP_1, sparse(femregion.ndof_p,femregion.ndof_p); ...
    sparse(femregion.ndof_p,femregion.ndof_p), A.MPrjP_1];

% Elastic dg matrix
V    = [A.V1,    A.V2;    A.V3,    A.V4];
IT_P = [A.IT1_P, A.IT2_P; A.IT3_P, A.IT4_P];
S_P  = [A.S1_P,  A.S2_P;  A.S3_P,  A.S4_P];

S_P_EP  = [A.S1_P_EP,  A.S2_P_EP;  A.S3_P_EP,  A.S4_P_EP];


% Poroelastic dg matriz
B   = [A.B1,   A.B2;   A.B3,   A.B4];
BT  = [A.BT1,  A.BT2;  A.BT3,  A.BT4];
S_B = [A.S1_B, A.S2_B; A.S3_B, A.S4_B];

BT_EP  = [A.BT1_EP, A.BT2_EP;  A.BT3_EP,  A.BT4_EP];
S_B_EP = [A.S1_B_EP, A.S2_B_EP; A.S3_B_EP, A.S4_B_EP];


B_beta   = [A.B1_beta, A.B2_beta; A.B3_beta, A.B4_beta];
BT_beta  = [A.BT1_beta, A.BT2_beta; A.BT3_beta, A.BT4_beta];
S_B_beta = [A.S1_B_beta, A.S2_B_beta; A.S3_B_beta, A.S4_B_beta];

BT_beta_EP  = [A.BT1_beta_EP,  A.BT2_beta_EP;  A.BT3_beta_EP,  A.BT4_beta_EP];
S_B_beta_EP = [A.S1_B_beta_EP, A.S2_B_beta_EP; A.S3_B_beta_EP, A.S4_B_beta_EP];


B_beta2   = [A.B1_beta2,   A.B2_beta2;   A.B3_beta2,   A.B4_beta2];
BT_beta2  = [A.BT1_beta2,  A.BT2_beta2;  A.BT3_beta2,  A.BT4_beta2];
S_B_beta2 = [A.S1_B_beta2, A.S2_B_beta2; A.S3_B_beta2, A.S4_B_beta2];

BT_beta2_EP  = [A.BT1_beta2_EP,  A.BT2_beta2_EP;  A.BT3_beta2_EP,  A.BT4_beta2_EP];
S_B_beta2_EP = [A.S1_B_beta2_EP, A.S2_B_beta2_EP; A.S3_B_beta2_EP, A.S4_B_beta2_EP];

% Imperfect pore matrix
D = [A.D1, A.D2; A.D3, A.D4];

% absorbing matrix
ABC_UU =  [A.ABC_uu_1, A.ABC_uu_2; A.ABC_uu_3, A.ABC_uu_4];
ABC_UW =  [A.ABC_uw_1, A.ABC_uw_2; A.ABC_uw_3, A.ABC_uw_4];
ABC_WU =  [A.ABC_wu_1, A.ABC_wu_2; A.ABC_wu_3, A.ABC_wu_4];
ABC_WW =  [A.ABC_ww_1, A.ABC_ww_2; A.ABC_ww_3, A.ABC_ww_4];

Matrices = ...
    struct( 'MPrjP_1', A.MPrjP_1,...
    'A_E', V - IT_P - transpose(IT_P) + S_P + S_P_EP,...
    'A_P',         B       - BT       - transpose(BT)       + S_B       - BT_EP - transpose(BT_EP) + S_B_EP,...
    'A_P_beta_pf', B_beta  - BT_beta  - transpose(BT_beta)  + S_B_beta  - BT_beta_EP            - (1-Data.delta)*transpose(BT_beta_EP)  + (1-Data.delta)*S_B_beta_EP,...
    'A_P_beta_fp', B_beta  - BT_beta  - transpose(BT_beta)  + S_B_beta  - (1-Data.delta)*BT_beta_EP  - transpose(BT_beta_EP)            + (1-Data.delta)*S_B_beta_EP,...
    'A_P_beta2',   B_beta2 - BT_beta2 - transpose(BT_beta2) + S_B_beta2 - (1-Data.delta)*BT_beta2_EP - (1-Data.delta)*transpose(BT_beta2_EP) + (1-Data.delta)^2*S_B_beta2_EP,...
    'MPrjP',MprjP,...
    'DGe', V + S_P,...
    'DGp',B + S_B ,...
    'DGp_beta',B_beta + S_B_beta ,...
    'M_P_rho', M_P_rho, ...
    'M_P_rhow', M_P_rhow, ...
    'M_P_rhof', M_P_rhof, ...
    'M_P_eta_kper', M_P_eta_kper, ...
    'D_vel', D_vel, ...
    'D_dis', D_dis, ...
    'D_P', D, ...
    'P1',A.P1,...
    'P2',A.P2,...
    'P1_beta',A.P1_beta,...
    'P2_beta',A.P2_beta,...
    'ABC_UU',ABC_UU,...
    'ABC_UW',ABC_UW,...
    'ABC_WU',ABC_WU,...
    'ABC_WW',ABC_WW);

fprintf('\n');

