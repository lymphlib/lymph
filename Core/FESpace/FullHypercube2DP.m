%> @file  FullHypercube2DP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief Evaluate the 2D modal basis and its derivatives hypercube at \f$(a,b)\f$ of order \f$(i,j)\f$
%> 
%> The function evaluates the 2D modal basis and its derivatives hypercube at \f$(a,b)\f$
%> of order \f$(i,j)\f$.
%>
%> 
%======================================================================
%> @section classFullHypercube2DP Class description
%======================================================================
%> @brief Evaluate the 2D modal basis and its derivatives hypercube at
%> \f$(a,b)\f$ of order \f$(i,j)\f$.
%>
%>
%> @param x     x-coordinates where evaluate the polynomials.
%> @param y     y-coordinates where evaluate the polynomials.
%> @param i     Orders for the construction of the polynomials in x-direction.
%> @param j     Orders for the construction of the polynomials in y-direction.
%>
%> @retval P    2D polynomial on the hypercube.
%> @retval dPdx    Derivative in x-direction of 2D polynomial on the hypercube.
%> @retval dPdy    Derivative in y-direction of 2D polynomial on the hypercube.
%======================================================================

function [P, dPdx, dPdy] = FullHypercube2DP(x, y, i, j)

    Lx  = LegendreP(x,i);
    Ly  = LegendreP(y,j);
    
    c = sqrt((2*i+1).*(2*j+1)/4);
    
    % psi
    P = c.*Lx.*Ly;
    
    if nargout == 3
    
        dLx = GradLegendreP(x, i);
        dLy = GradLegendreP(y, j);
        
        % x-derivative
        dPdx = c.*dLx.*Ly;
        
        % y-derivative
        dPdy = c.*Lx.*dLy;
    
    end

end
