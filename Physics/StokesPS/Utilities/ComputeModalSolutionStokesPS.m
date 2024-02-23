%> @file  ComputeModalSolutionStokesPS.m
%> @author Ilario Mazzieri
%> @date 16 February 2024
%> @brief Compute modal expansions for pseudo stress sigma
%>
%==========================================================================
%> @section ComputeModalSolutionStokesPS Class description
%> @brief  Compute Modal Solution for the unsteady Stokes problem
%
%> @param Data        Struct with problem's data
%> @param femregion   Finite Element struct (see CreateDOF.m)
%
%> @retval sigma       Vector with modal expansions
%>
%==========================================================================

function [sigma] = ComputeModalSolutionStokesPS(Data,femregion)
%% Quadrature values
[~, ~, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization
sigma11 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
sigma12 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
sigma21 = spalloc(femregion.ndof_f,1,femregion.ndof_f);
sigma22 = spalloc(femregion.ndof_f,1,femregion.ndof_f);


%% Loop over the elements
% Visualization of computational progress
index_shift=0;
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie=1:femregion.nel

    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);

    % Element tag
    tag_ie = femregion.tag(ie);

    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    index_element = index_shift + (1:1:femregion.nedges(ie))';
    index_shift   = index_element(end);

    % Extraction of element geometrical information
    coords_ie          = femregion.coords_element{ie};

    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

    % Check if the element is fluid
    if tag_ie == 'F'

        index_f = index - femregion.ndof_p;
        % Inizialize local silution vector
        s11_loc = 0;
        s12_loc = 0;
        s21_loc = 0;
        s22_loc = 0;

        for iTria = 1:size(Tria,1)

            % Construction of Jacobian and quadrature nodes
            [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);

            % Scaled weights
            dx = det(BJ) * w_2D;

            % Construction of the basis functions
            [phiq, ~, ~] = Evalshape2D(femregion, ie, qNodes_2D);

            % Compute exact solution
            s11ex_loc = Data.sigma_ex{1}(xq,yq);
            s12ex_loc = Data.sigma_ex{2}(xq,yq);
            s21ex_loc = Data.sigma_ex{3}(xq,yq);
            s22ex_loc = Data.sigma_ex{4}(xq,yq);

            s11_loc = s11_loc + (dx.*phiq)'*s11ex_loc;
            s12_loc = s12_loc + (dx.*phiq)'*s12ex_loc;
            s21_loc = s21_loc + (dx.*phiq)'*s21ex_loc;
            s22_loc = s22_loc + (dx.*phiq)'*s22ex_loc;

        end

        sigma11(index_f) = s11_loc;
        sigma12(index_f) = s12_loc;
        sigma21(index_f) = s21_loc;
        sigma22(index_f) = s22_loc;

    end

end

sigma = [sigma11; sigma12; sigma21; sigma22];

fprintf('\n');

