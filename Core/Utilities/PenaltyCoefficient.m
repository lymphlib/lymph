%> @file  PenaltyCoefficient.m
%> @author Stefano Bonetti, Mattia Corti, Ivan Fumagalli, Ilario Mazzieri,
%> Caterina Leimer Saglio
%> @date 10 April 2026
%> @brief Construction of the penalty coefficients for all the faces of an
%> element.
%> 
%> The function constructs the penalty coefficients for all the faces of an
%> element, only considering the geometrical component of the penalty. This
%> function does NOT provide any multiplication by physical parameters in
%> the PDE. 
%> 
%======================================================================
%> @section classPenaltyCoefficient Class description
%======================================================================
%> @brief  Construction of the penalty coefficients for all the faces of an
%> element.
%> @param femregion      Structure containing the information about FEM 
%> discretization.
%> @param Data           Structure containing the information about the
%> discretization and the physical parameters of the problem.
%> @param ie             ID of the evaluated mesh element.
%> @param neighedges_ie  ID of the edges of the mesh element, with respect
%> to the enumeration of the neighboring element.
%> @param neigh_ie       ID of the neighboring elements.
%> @param meshsize       Diameters of the faces of the element (length).
%>
%> @retval penalty     Struct containing the arrays with the geometrical part of the penalty
%> coefficient for all the faces of an element (on all the boundary faces the 
%> penalty coefficient is defined as in the Dirichlet case). The
%> implemented policies are:
%>   \c max \cite cangiani2017,
%>   \c min \cite antonietti2022stability,
%>   or harmonic average - \c harm - \cite dryja2007bddc between the elements.
%======================================================================

function [penalty] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize)
    
    %% Creation of the mask of internal edges
    neigh_app = neigh_ie(neigh_ie>0);
    neighedges_app = neighedges_ie(neigh_ie>0);
    mask = diag(neigh_ie>0);
    mask = mask(any(mask==1),:);

    %% Construction of penalty part dependent on the femregion degree
    deg = femregion.degree.*(femregion.degree>0)+(femregion.degree==0);

    %% Construction of inverse trace constant
    Cinv = femregion.area(ie)./femregion.max_kb{ie};

    %% External Cinv computation
    Cinv_ext = zeros(size(neigh_app))';
    
    for ii = 1:length(neigh_app)
        Cinv_ext(ii) = femregion.area(neigh_app(ii))./femregion.max_kb{neigh_app(ii)}(neighedges_app(ii));
    end
    
    %% Penalty max construction
    p_max   = (Data.penalty_coeff.* deg(ie)^2 .* Cinv .* meshsize/femregion.area(ie))';
    p_max_n = (Data.penalty_coeff.* ((deg(neigh_app).^2.*Cinv_ext./femregion.area(neigh_app))'*mask)' .* meshsize)';
    penalty.max = max(p_max, p_max_n);

    %% Penalty min construction
    p_min   = (Data.penalty_coeff./deg(ie) .* femregion.area(ie)./meshsize)';
    p_min_n = (Data.penalty_coeff .* ((femregion.area(neigh_app)./deg(neigh_app))'*mask)'./meshsize)';
    penalty.min = (neigh_ie>0).*min(p_min,p_min_n) + (neigh_ie<0).*p_min;

    %% Penalty harmonic construction
    p_harm   = (Data.penalty_coeff.* deg(ie)^2 .* Cinv .* meshsize/femregion.area(ie))';
    p_harm_n = (Data.penalty_coeff.* ((deg(neigh_app).^2.*Cinv_ext./femregion.area(neigh_app))'*mask)' .* meshsize)';
    penalty.harm = (neigh_ie>0).*mean([p_harm; p_harm_n]) + (neigh_ie<0).*p_harm;

end

