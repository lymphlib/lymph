%> @file  JacobiP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 16 September 2025 
%> @brief The function evaluates the Jacobi polynomials
%> 
%> The function evaluates the Jacobi polynomial with parameters \f$\alpha\f$ and \f$\beta\f$ 
%> of order \f$N\f$ at point \f$x\f$, namely \f$P_N^{(\alpha,\beta)}(x)\f$. The polynomials are
%> hard-coded until the degree 1. Then we compute them by means of the recursive formula \cite jacobiformulas.
%>

%> Note: They are normalized to be orthonormal.
%> 
%======================================================================
%> @section classJacobiP Class description
%======================================================================
%> @brief The function evaluates the Jacobi polynomials.
%>
%> @param x     Vector of points \f$x_i\f$ in which the Jacobi polynomials are evaluated.
%> @param alpha First numerical parameter of Jacobi polynomials.
%> @param beta  Second numerical parameter of Jacobi polynomials.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval P   Matrix containing at the position \f$P_{ij}\f$ a Jacobi polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$, computed with parameters \f$(\alpha,\beta)\f$.
%===================================================s===================

function [P] = JacobiP(x, alpha, beta, N)

    % Turn points into row if needed
    x = x(:); 

    % Extract maximum degree of the polinomial
    Nmax = max(N);

    % Construction of initial structure for polynomial construction
    Ptable  = zeros(length(x), Nmax + 1);

    % Polynomial constants
    gamma0 = (2^(alpha+beta+1)/(alpha+beta+1)) * ... 
        gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+1);

    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;

    % Construction of the basis for constant function
    Ptable(:,1) = ones(length(x),1)/sqrt(gamma0);
    
    % Construction of the basis for the linear function
    Ptable(:,2) = (((alpha+beta+2)*x/2 + (alpha-beta)/2)/sqrt(gamma1));
    
    %% Recursive formula for computation of polynomials
    c_oold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));
    
    for ii = 1:Nmax-1

        % Compute the constants of the recursive formula
        c_app = 2*ii+alpha+beta;
        c_new = 2/(c_app+2)*sqrt((ii+1)*(ii+1+alpha+beta)* ...
            (ii+1+alpha)*(ii+1+beta)/(c_app+1)/(c_app+3));
        c_old = - (alpha^2-beta^2)/(c_app*(c_app+2));

        Ptable(:, ii+2)  = 1/c_new*((x-c_old).*Ptable(:, ii+1) - c_oold.*Ptable(:, ii));

        c_oold = c_new;
    end

    % Output construction
    P  = Ptable(:, N + 1);
    
  end
