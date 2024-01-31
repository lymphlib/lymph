%> @file  PenaltyCoefficient.m
%> @author Stefano Bonetti, Mattia Corti, Ivan Fumagalli, Ilario Mazzieri
%> @date 29 March 2023 
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
%> @param penalty_policy [max,min,harm] (optional: default=max) Choice for
%>                       penalty coefficient, that can be based on either
%>                       max \cite cangiani2017,
%>                       min \cite antonietti2022stability,
%>                       or harmonic average \cite dryja2007bddc
%>                       between current and neighboring element.
%>
%> @retval penalty     Array containing the geometrical part of the penalty
%> coefficient for all the faces of an element.
%======================================================================

function [penalty] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize, penalty_policy)

    if ~exist('penalty_policy','var')
        penalty_policy = 'max';
    end

    % Treatment of the degree 0 case
    if femregion.degree == 0 
        penalty_coeff = Data.penalty_coeff;
    else
        if strcmp(penalty_policy,'max')
            penalty_coeff = Data.penalty_coeff.*(femregion.degree^2);
        elseif strcmp(penalty_policy,'min')
            penalty_coeff = Data.penalty_coeff./femregion.degree;
        elseif strcmp(penalty_policy,'harm')
            %> @TODO check dependency on degree in the harm case
            penalty_coeff = Data.penalty_coeff./femregion.degree;
        else
            error(strcat('Unknown penalty policy: ', penalty_policy))
        end
    end

    % Construction of inverse trace constant
    Cinv = femregion.area(ie)./femregion.max_kb{ie};
    
    % Dirichlet boundary face
    penalty_dir = [];
    if strcmp(penalty_policy,'max')
        penalty_dir = penalty_coeff * Cinv .* (meshsize/femregion.area(ie));
    elseif strcmp(penalty_policy,'min')
        penalty_dir = penalty_coeff * (femregion.area(ie)/meshsize);
    elseif strcmp(penalty_policy,'harm')
        penalty_dir = penalty_coeff * (meshsize/femregion.area(ie));
    else
        error(strcat('Unknown penalty policy: ', penalty_policy))
    end

    % Internal face: value of the penalty for the element itself
    p_ie = zeros(size(neigh_ie))';

    % Internal face: value of the penalty for the neighboring element
    Cinv_ext = zeros(size(neigh_ie));
    p_ie_n = zeros(size(neigh_ie))';
    
    for ii = 1:length(neigh_ie)
        % Computation for the neighbor elements
        if neigh_ie(ii) > 0
            if strcmp(penalty_policy,'max')
                Cinv_ext(ii) = femregion.area(neigh_ie(ii))./femregion.max_kb{neigh_ie(ii)}(neighedges_ie(ii));
                p_ie(ii) = penalty_coeff * Cinv(ii) .* (meshsize(ii)/femregion.area(ie));
                p_ie_n(ii) = penalty_coeff * Cinv_ext(ii) * (meshsize(ii)/femregion.area(neigh_ie(ii)));
            elseif strcmp(penalty_policy,'min')
                p_ie(ii) = penalty_coeff * femregion.area(ie)/meshsize(ii);
                p_ie_n(ii) = penalty_coeff * femregion.area(neigh_ie(ii))/meshsize(ii);
            elseif strcmp(penalty_policy,'harm')
                p_ie(ii) = penalty_coeff * (meshsize(ii)/femregion.area(ie));
                p_ie_n(ii) = penalty_coeff * femregion.area(neigh_ie(ii))/meshsize(ii);
            else
                error(strcat('Unknown penalty policy: ', penalty_policy))
            end
        end

    end
    
    % Computation of the penalty for all the faces of the element
    if strcmp(penalty_policy,'max')
        penalty = max([p_ie p_ie_n],[],2);
    elseif strcmp(penalty_policy,'min')
        penalty = penalty_coeff.*min([p_ie p_ie_n],[],2);
    elseif strcmp(penalty_policy,'harm')
        % Note that in \cite dryja2007bddc , l_ij=2 is there because the sum in (6) is
        % over the faces of an element, thus containing only half of the contribution.
        % On the other hand, here we do not need l_ij.
        penalty = penalty_coeff.*(p_ie + p_ie_n)./(2*p_ie.*p_ie_n);
    else
        error(strcat('Unknown penalty policy: ', penalty_policy))
    end
    %> @TODO check the cases neigh_ie<-1 and introduce comments to explain the choice
    % penalty = (neigh_ie > 0)'.*max([p_ie p_ie_n],[],2) + (neigh_ie == -1)'.*penalty_dir;
    penalty = (neigh_ie > 0)'.*penalty + (neigh_ie == -1)'.*penalty_dir;

end

