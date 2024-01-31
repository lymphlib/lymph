%> @file  Simplex2DP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief Evaluate 2D orthonormal polynomial on simplex at \f$(a,b)\f$ of order \f$(i,j)\f$
%> 
%> The function evaluates 2D orthonormal polynomial on simplex at 
%> \f$(a,b)\f$ of order \f$(i,j)\f$.
%>
%> 
%======================================================================
%> @section classSimplex2DP Class description
%======================================================================
%> @brief Evaluate 2D orthonormal polynomial on simplex at \f$(a,b)\f$ of
%> order \f$(i,j)\f$.
%>
%> @param a     x-coordinates where evaluate the polynomials.
%> @param b     y-coordinates where evaluate the polynomials.
%> @param i     Orders for the construction of the polynomials in x-direction.
%> @param j     Orders for the construction of the polynomials in y-direction.
%>
%> @retval P    2D polynomial on the simplex.
%======================================================================

function [P] = Simplex2DP(a, b, i, j)

    % Computation of 1D Jacobi polynomials
    Lx = JacobiP(a, 0, 0, i);
    Ly = JacobiP(b, 2*i+1, 0, j);

    % Computation of 2D Jacobi polynomials
    P = sqrt(2.0)*Lx.*Ly.*(1-b).^i;

return;
