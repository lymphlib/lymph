%> @file  GradLegendreP.m
%> @author Mattia Corti, Paola F. Antonietti, Matteo Caldana, Caterina Leimer Saglio
%> @date 28 October 2025 
%> @brief The function evaluates Legendre polynomials and their first
%> derivatives.
%> 
%> The function evaluates the Legendre polynomials and its gradient. 
%> We compute them by means of the recursive formula \cite handbookformulas :
%> \f[ P_{i+1}(x) = \frac{1}{i+1}\{(2i+1)xP_i(x)-iP_{i-1}(x)\} \f]
%> \f[ P'_{i+1}(x) = x P'_{i}(x) + (i+1)P_{i}(x) \f]
%>
%======================================================================
%> @section classGradLegendreP Class description
%======================================================================
%> @brief The function evaluates Legendre polynomials and their first derivatives.
%>
%> @param x     Vector of points \f$x_i\f$ in which the Legendre polynomials are evaluated.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval P    Matrix containing at the position \f$P_{ij}\f$ a Legendre polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$.
%> @retval dP   Matrix containing at the position \f$P_{ij}'\f$ a gradient of a Legendre polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$.
%======================================================================

function [P, dP] = GradLegendreP(x,N)
    
    % Turn points into column if needed
    x = x(:); 

    % Extract maximum degree of the polinomial
    Nmax = max(N);

    % Construction of initial structure for polynomial construction
    Ptable   = ones(length(x), Nmax + 1);
    dPtable  = zeros(length(x), Nmax + 1);
    
    Ptable(:, 2)   = x;
    dPtable(:, 2)  = ones(length(x), 1);
    
    %% Recursive formula for computation of polynomials
    for i = 1:Nmax-1
        Ptable(:, i+2)   = ((2*i + 1)*x.*Ptable(:, i+1) - i*Ptable(:, i))/(i + 1);
        dPtable(:, i+2)  = (i + 1) * Ptable(:, i+1) + x .* dPtable(:, i+1);
    end

    % Output construction
    P   = Ptable(:, N + 1);
    dP  = dPtable(:, N + 1);
    
end
