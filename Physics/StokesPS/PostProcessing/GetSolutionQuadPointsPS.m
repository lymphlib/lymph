%> @file  GetSolutionQuadPointsPS.m
%> @author Ilario Mazzieri
%> @date 16 February 2024
%> @brief  Compute solution at quadrature nodes
%>
%==========================================================================
%> @section classGetSolutionQuadPointsPS Class description
%==========================================================================
%> @brief  Compute solution at quadrature nodes
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Solutions  Struct with solution vectors
%> @param time       time instant
%
%> @retval Displacement   Vector with solution at quad nodes
%> @retval Pressure       Vector with solution pressure at quad nodes
%>
%==========================================================================

function [Displacement, Pressure] = GetSolutionQuadPointsPS(Data, femregion, Solutions, time)

%% Quadrature values
[~, ~, ref_qNodes_2D, ~] = Quadrature(femregion.nqn);

%% Setup
% G = [ x | y | sigma_h ... | sigma_ex ... ];

% Pseudo-stree
Sigma  = zeros(femregion.nel_f*femregion.nqn.^2,10);
% Acceleration
divSigma  = zeros(femregion.nel_f*femregion.nqn.^2,4);
% Pressure
Pf   = zeros(femregion.nel_f*femregion.nqn.^2,4);

% shift index
kindf = 0;


%% Loop over the elements
% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie = 1:femregion.nel % loop over elements

    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);

    % Element id
    id_ie  = femregion.id(ie);

    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

    % Extraction of element geometrical information
    coords_ie = femregion.coords_element{ie};

    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

    % Fluid element
    if femregion.tag(ie) == 'F'

        index_f = index;
        if time == 0
            % Local solution
            s11_loc = Solutions.sigma_h(index_f);
            s12_loc = Solutions.sigma_h(index_f + 1*femregion.ndof_f);
            s21_loc = Solutions.sigma_h(index_f + 2*femregion.ndof_f);
            s22_loc = Solutions.sigma_h(index_f + 3*femregion.ndof_f);

        else
            % Local solution
            s11_loc = Solutions{end}.sigma_h(index_f);
            s12_loc = Solutions{end}.sigma_h(index_f + 1*femregion.ndof_f);
            s21_loc = Solutions{end}.sigma_h(index_f + 2*femregion.ndof_f);
            s22_loc = Solutions{end}.sigma_h(index_f + 3*femregion.ndof_f);
        end


        for iTria = 1:size(Tria,1)

            % Construction of Jacobian and quadrature nodes
            [~, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

            xq  = qNodes_2D(:,1);
            yq  = qNodes_2D(:,2);
            lqn = length(xq);

            % Construction of the basis functions
            [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

            % Exact and approximated solutions at quadrature points
            s11ex_loc = Data.sigma_ex{1}(xq,yq)*Data.sigma_t_ex{1}(time);
            s12ex_loc = Data.sigma_ex{2}(xq,yq)*Data.sigma_t_ex{1}(time);
            s21ex_loc = Data.sigma_ex{3}(xq,yq)*Data.sigma_t_ex{1}(time);
            s22ex_loc = Data.sigma_ex{4}(xq,yq)*Data.sigma_t_ex{1}(time);

            s11h_loc  = phiq*s11_loc;
            s12h_loc  = phiq*s12_loc;
            s21h_loc  = phiq*s21_loc;
            s22h_loc  = phiq*s22_loc;

            divs1h_loc = gradqx*s11_loc + gradqy*s12_loc;
            divs2h_loc = gradqx*s21_loc + gradqy*s22_loc;


            % Fill the output Phi and Pa
            Sigma(kindf + 1 : kindf + lqn,:) = [xq, yq, s11h_loc, s12h_loc, s21h_loc, s22h_loc, ...
                s11ex_loc, s12ex_loc, s21ex_loc, s22ex_loc];

            % a = div(sigma)
            divSigma(kindf + 1 : kindf + lqn,:) = [xq, yq, divs1h_loc, divs2h_loc];

            % p = -0.5 tr(sigma)
            Pf(kindf + 1 : kindf + lqn,:) = [xq, yq, -0.5*(s11h_loc + s22h_loc), -0.5*(s11ex_loc + s22ex_loc)];


            kindf = kindf + lqn;

        end

    end

end

%% Output
Displacement.Sigma = Sigma;
Displacement.SigmaTag = '\sigma';
Displacement.divSigma = divSigma;
Displacement.divSigmaTag = 'div(\sigma)';

Pressure.Pf = Pf;


