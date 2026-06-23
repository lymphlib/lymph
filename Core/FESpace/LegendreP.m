%> @file  LegendreP.m
%> @author Mattia Corti, Paola F. Antonietti, Matteo Caldana
%> @date 28 October 2025 
%> @brief The function evaluates Legendre polynomials.
%> 
%> The function evaluates the Legendre polynomials. We compute them by 
%> means of the recursive formula \cite handbookformulas :
%> \f[ P_{i+1}(x) = \frac{1}{i+1}\{(2i+1)xP_i(x)-iP_{i-1}(x)\} \f]
%>
%======================================================================
%> @section classLegendreP Class description
%======================================================================
%> @brief The function evaluates Legendre polynomials.
%>
%> @param x     Vector of points \f$x_i\f$ in which the Legendre polynomials are evaluated.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval P    Matrix containing at the position \f$P_{ij}\f$ a Legendre polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$.
%======================================================================

function [P] = LegendreP(x, N)
    
    % Turn points into column if needed
    x = x(:); 

    % Extract maximum degree of the polinomial
    Nmax = max(N);

    % Construction of initial structure for polynomial construction
    Ptable  = ones(length(x), Nmax + 1);
    
    Ptable(:, 2)  = x;

    %% Recursive formula for computation of polynomials
    for i = 1:Nmax-1
        Ptable(:, i+2)  = ((2*i+1)*x.*Ptable(:, i+1) - i*Ptable(:, i))/(i + 1);
    end

    % Output construction
    P  = Ptable(:, N + 1);
end
