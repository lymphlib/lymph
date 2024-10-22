%> @file  Evalshape2DCoeff.m
%> @author Mattia Corti
%> @date 13 November 2023 
%> @brief Construction of the coefficients of the basis functions.
%> 
%> The function constructs the coefficients of the polynomial basis functions \f$\varphi_j(x,y)=\mathcal{L}_i(x)\mathcal{L}_k(y)\f$.
%> 
%======================================================================
%> @section classEvalshape2DCoeff Class description
%======================================================================
%> @brief  Construction of the the structures for the basis functions.
%>
%> @param N             Maximum degree of the basis functions. 
%>
%> @retval Lx          Coefficients of the Legendre polynomials in x-direction.
%> @retval Ly          Coefficients of the Legendre polynomials in y-direction.
%> @retval dLx         Coefficients of the gradient of Legendre polynomials in x-direction.
%> @retval dLy         Coefficients of the gradient of Legendre polynomials in y-direction.
%======================================================================

function [Lx, Ly, dLx, dLy] = Evalshape2DCoeff(N)

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
    
    % Construction of the coefficients of the Legendre polynomials derivatives in 1D
    dLx = Lx.*(0:size(Lx,1)-1)';
    dLx = [dLx(2:end,:); dLx(1,:)];

    dLy = Ly.*(0:size(Ly,1)-1)';
    dLy = [dLy(2:end,:); dLy(1,:)];

end
