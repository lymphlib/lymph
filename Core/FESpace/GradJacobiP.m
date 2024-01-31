%> @file  GradJacobiP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief The function evaluates the gradient of Jacobi polynomials
%> 
%> The function evaluates the gradient of Jacobi polynomial with parameters \f$\alpha\f$ and \f$\beta\f$ 
%> of order \f$N\f$ at point \f$x\f$, namely \f$P_N^{(\alpha,\beta)}(x)\f$. The polynomials are
%> computed by using the relation:
%> \f[ {P'}_{N}^{(\alpha,\beta)}(x) = \sqrt{N(N+\alpha+\beta+1)} P_{N-1}^{(\alpha+1,\beta+1)}(x) \f]
%>

%> Note: They are normalized to be orthonormal.
%> 
%======================================================================
%> @section classGradJacobiP Class description
%======================================================================
%> @brief The function evaluates the gradient of Jacobi polynomials.
%>
%> @param x     Vector of points \f$x_i\f$ in which the gradient of Jacobi polynomials are evaluated.
%> @param alpha First numerical parameter of gradient of Jacobi polynomials.
%> @param beta  Second numerical parameter of gradient of Jacobi polynomials.
%> @param N     Vector of gradient of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval dP   Matrix containing at the position \f$P_{ij}'\f$ a gradient of a Jacobi polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$, computed with parameters \f$(\alpha,\beta)\f$.
%===================================================s===================

function [dP] = GradJacobiP(x, alpha, beta, N)

    % Computation of the associated Jacobi polynomials
    dP = sqrt(N*(N+alpha+beta+1))*JacobiP(x, alpha+1, beta+1, N-1);

    % Derivative of the constant polynomial is equal to 0
    if any(N == 0)
        dP(N == 0) = 0.0;
    end

end
