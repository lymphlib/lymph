%> @file  Evalshape2DCoeff.m
%> @author Mattia Corti
%> @date 13 November 2023 
%> @brief Construction of the coefficients of the basis functions for a neighboring element.
%> 
%> The function constructs the coefficients of the polynomial basis functions 
%> \f$\varphi_j(x,y)=\mathcal{L}_i(x)\mathcal{L}_k(y)\f$ for a neighboring 
%> element. For the details of formulas, we refer to \cite pennesiquadraturefree
%> 
%======================================================================
%> @section classDecr Class description
%======================================================================
%> @brief  Construction of the the structures for the basis functions for a neighboring element.
%>
%> @param N            Maximum degree of the basis functions. 
%> @param BJ_tilde     Transformation matrix BJ_tilde.
%> @param tr_tilde     Translation vector tr_tilde.
%> @param coeff_bin    Matrix of binomial coefficients. 
%>
%> @retval Lx_neigh    Coefficients of the Legendre polynomials in x-direction.
%> @retval Ly_neigh    Coefficients of the Legendre polynomials in y-direction.
%======================================================================


function [Lx_neigh, Ly_neigh] = Evalshape2DCoeffNeigh(N, BJ_tilde, tr_tilde, coeff_bin)

    % First indexes vector
    idx1 = reshape(tril((1:N+1).*ones(N+1,N+1)),[1 (N+1)^2]);
    idx1(idx1==0) = [];
    idx1 = idx1 - 1;

    % Second indexes vector 
    idx2 = reshape(flip(triu((N+1:-1:1).*ones(N+1,N+1)),2)',[1 (N+1)^2]);
    idx2(idx2==0) = [];
    idx2 = idx2 - 1;

    % Construction of the coefficients of the Legendre polynomials in 1D
    Lx  = LegendrePCoeff(idx1).*sqrt((2*idx1+1)/2);
    Ly  = LegendrePCoeff(idx2).*sqrt((2*idx2+1)/2);

    Lx_neigh = zeros(size(Lx));
    Ly_neigh = zeros(size(Ly));

    % Construction of the neighbours terms
    for jj = 0:size(Lx,1)-1
        for ii = 0:jj
            Lx_neigh(ii+1,:) = Lx_neigh(ii+1,:) + Lx(jj+1,:)*(coeff_bin(ii+1,jj+1)'.*BJ_tilde(1,1).^ii.*tr_tilde(1).^(jj-ii));
            Ly_neigh(ii+1,:) = Ly_neigh(ii+1,:) + Ly(jj+1,:)*(coeff_bin(ii+1,jj+1)'.*BJ_tilde(2,2).^ii.*tr_tilde(2).^(jj-ii));
       end
    end

end
