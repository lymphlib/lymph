%> @file  CoeffBin.m
%> @author Mattia Corti
%> @date 11 November 2023 
%> @brief Function that construct a matrix of binomial coefficients.
%> 
%> The function constructs the matrix of binomial coefficients, in the
%> sense that in position (i,j) we found the binomial coefficient j over i.
%>
%======================================================================
%> @section classDecr Class description
%======================================================================
%> @brief Function that construct a matrix of binomial coefficients.
%>
%> @param deg           Maximum degree of binomial coefficient to compute,
%> 
%> @retval coeff_bin    Matrix of binomial coefficients.
%======================================================================

function [coeff_bin] = CoeffBin(deg)

    % Preallocation of binomial coefficients matrix
    coeff_bin = zeros(deg+1,deg+1);

    % Cycle over the first term of the binomial coefficient
    for jj = 0:deg

        % Cycle over the second term of the binomial coefficient (ii<=jj)
        for ii = 0:jj
            coeff_bin(ii+1,jj+1) = nchoosek(jj,ii);
        end
    
    end

end