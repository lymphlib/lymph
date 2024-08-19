%> @file  ComputeErrorInterfacePA.m
%> @author Ilario Mazzieri
%> @date 28 June 2024
%> @brief  Compute errors for the poro-elasto-acoustic problem (only
%interface between poro and acoustic material)
%>
%==========================================================================
%> @section classComputeErrorInterfacePA Class description
%> @brief  Compute errors for waves's problem
%
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param w_h        Struct with solution vectors
%> @param time       Parameter time instant
%> @param tau        Parameter in [0,1] to define interface poro condition
%
%> @retval Ew_interface     Computed error at the interface
%>
%==========================================================================
function [Ew_interface] = ComputeErrorInterfacePA(Data,femregion,neighbor,w_h,time,tau)

% initialization
Ew_interface = 0;

%% Quadrature values
[ref_qNodes_1D, w_1D, ~, ~] = Quadrature(femregion.nqn);

for ie = 1 : femregion.nel

    % Id and tag selection for elements
    tag_ie = femregion.tag(ie);

    % Selection of the matrix positions associated to element ie
    index         = (ie-1)*femregion.nln*ones(femregion.nln,1) + (1:femregion.nln)';

    % Extraction of neighbor element and their edges
    neigh_ie      = neighbor.neigh{ie};

    % Extraction of element geometrical information
    coords_ie          = femregion.coords_element{ie};
    [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);

    if tag_ie == 'P'

        for iedg = 1: neighbor.nedges(ie) % loop over faces

            if (neigh_ie(iedg) >0 && femregion.tag(neigh_ie(iedg)) == 'A')

                % Extraction filtration displacement
                local_wh1 = w_h(index);
                local_wh2 = w_h(index+femregion.ndof_p);

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

                % Construction of the basis functions
                phiedgeq = Evalshape2D(femregion, ie, qNodes_1D);

                local_exact_w1 = Data.wp_ex{1}(xq,yq).*Data.wp_t_ex{1}(time)';
                local_exact_w2 = Data.wp_ex{2}(xq,yq).*Data.wp_t_ex{1}(time)';

                local_aprox_w1 = local_wh1*phiedgeq;
                local_aprox_w2 = local_wh2*phiedgeq;


                Ew_interface = Ew_interface + (1-tau)/tau * ds *  ...
                    ((local_aprox_w1 - local_exact_w1) * nx ...
                    + (local_aprox_w2 - local_exact_w2) * ny).^2;

            end
        end
    end
end

