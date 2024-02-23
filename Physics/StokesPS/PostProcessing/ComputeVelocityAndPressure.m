%> @file  ComputeVelocityAndPressure.m
%> @author Ilario Mazzieri
%> @date 16 Febraury 2024
%> @brief  Compute velocity and pressure field form pseudo-stress
%>
%==========================================================================
%> @section classComputeVelocityAndPressure Class description
%> @brief  Compute velocity and pressure field form pseudo-stress
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Solutions  Struct with solution vector
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function [Velf,Pf] = ComputeVelocityAndPressure(Data, femregion, Solutions)

%% Quadrature values
[~, ~, ref_qNodes_2D, ~] = Quadrature(femregion.nqn);

%% Setup
% Pressure
Pf   = zeros(femregion.nel_f*femregion.nqn.^2,3);

% Velocity
divSigma  = zeros(femregion.nel_f*femregion.nqn.^2,4,size(Solutions,2)+1);


for i = 1 : size(Solutions,2)

    % shift index
    kindf = 0;

    %% Loop over the elements

    for ie = 1:femregion.nel % loop over elements

        % Element id
        %id_ie  = femregion.id(ie);

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

            index_f = index - femregion.ndof_p;

            % Local solution
            s11_loc = Solutions{i}.sigma_h(index_f);
            s12_loc = Solutions{i}.sigma_h(index_f + 1*femregion.ndof_f);
            s21_loc = Solutions{i}.sigma_h(index_f + 2*femregion.ndof_f);
            s22_loc = Solutions{i}.sigma_h(index_f + 3*femregion.ndof_f);


            for iTria = 1:size(Tria,1)

                % Construction of Jacobian and quadrature nodes
                [~, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

                xq  = qNodes_2D(:,1);
                yq  = qNodes_2D(:,2);
                lqn = length(xq);

                % Construction of the basis functions
                [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

                % Approximated solutions at quadrature points
                s11h_loc  = phiq*s11_loc;
                %s12h_loc  = phiq*s12_loc;
                %s21h_loc  = phiq*s21_loc;
                s22h_loc  = phiq*s22_loc;

                % Divergence at quadrature nodes
                divs1h_loc = gradqx*s11_loc + gradqy*s12_loc;
                divs2h_loc = gradqx*s21_loc + gradqy*s22_loc;


                % a = div(sigma)
                divSigma(kindf + 1 : kindf + lqn,:,i+1) = [xq, yq, divs1h_loc, divs2h_loc];

                % p = -0.5 tr(sigma)
                Pf(kindf + 1 : kindf + lqn,:) = [xq, yq, -0.5*(s11h_loc + s22h_loc)];


                kindf = kindf + lqn;

            end

        end

    end

end

%%
Velf = zeros(size(divSigma(:,:,1),1),size(divSigma(:,:,1),2));
Velf(:,1:2) = divSigma(:,1:2,2);
xq = divSigma(:,1,2);
yq = divSigma(:,2,2);
for i = 1 : size(Solutions,2)
    Velf(:,3) = Velf(:,3) + (divSigma(:,3,i) + divSigma(:,3,i+1))*Data.dt*0.5 ...
        + (Data.source_vel{1}(xq,yq)*Data.source_vel_t{1}(i*Data.dt) + Data.source_vel{1}(xq,yq)*Data.source_vel_t{1}((i+1)*Data.dt))*Data.dt*0.5 ...
        + (Data.source_vel_d{1}(xq,yq)*Data.source_vel_d_t{1}(i*Data.dt) + Data.source_vel_d{1}(xq,yq)*Data.source_vel_d_t{1}((i+1)*Data.dt))*Data.dt*0.5;
    Velf(:,4) = Velf(:,4) + (divSigma(:,4,i) + divSigma(:,4,i+1))*Data.dt*0.5 ...
        + (Data.source_vel{2}(xq,yq)*Data.source_vel_t{1}(i*Data.dt) + Data.source_vel{2}(xq,yq)*Data.source_vel_t{1}((i+1)*Data.dt))*Data.dt*0.5 ...
        + (Data.source_vel_d{2}(xq,yq)*Data.source_vel_d_t{1}(i*Data.dt) + Data.source_vel_d{2}(xq,yq)*Data.source_vel_d_t{1}((i+1)*Data.dt))*Data.dt*0.5;

end

%% Output
Velf(:,3) = Velf(:,3) + Data.vel0{1}(xq,yq);
Velf(:,4) = Velf(:,4) + Data.vel0{2}(xq,yq);
