%> @file  FullHypercube2DP.m
%> @author Mattia Corti, Paola F. Antonietti, Matteo Caldana, Caterina Leimer Saglio
%> @date 28 October 2025
%> @brief Evaluate the 2D modal basis and its derivatives hypercube at \f$(x,y)\f$ of order \f$(i,j)\f$.
%> The Legendre polynomial basis is used over a tensor-product domain.
%>
%> 
%======================================================================
%> @section classFullHypercube2DP Class description
%======================================================================
%> @brief Evaluate the 2D modal basis and its derivatives hypercube at
%> \f$(x,y)\f$ of order \f$(i,j)\f$.
%> The Legendre polynomial basis is used over a tensor-product domain.
%>
%>
%> @param x     x-coordinates where evaluate the polynomials.
%> @param y     y-coordinates where evaluate the polynomials.
%> @param i     Orders for the construction of the polynomials in x-direction.
%> @param j     Orders for the construction of the polynomials in y-direction.
%>
%> @retval P       2D polynomial on the hypercube.
%> @retval dPdx    Derivative in x-direction of 2D polynomial on the hypercube.
%> @retval dPdy    Derivative in y-direction of 2D polynomial on the hypercube.
%> @retval d2Pdx2  Second derivative in x-direction of 2D polynomial on the hypercube.
%> @retval d2Pdy2  Second derivative in y-direction of 2D polynomial on the hypercube.
%> @retval d2Pdxy  Second mixed derivative of 2D polynomial on the hypercube.
%======================================================================

function [P, dPdx, dPdy, d2Pdx2, d2Pdy2, d2Pdxy] = FullHypercube2DP(x, y, i, j)

    c = sqrt((2*i+1).*(2*j+1)/4);
    
    switch nargout
        case 6

            [Lx, dLx, d2Lx] = LapLegendreP(x, i);
            [Ly, dLy, d2Ly] = LapLegendreP(y, j);

            % x-derivatives
            d2Pdx2 = c.*d2Lx.*Ly;
            dPdx   = c.*dLx.*Ly;

            % y-derivatives
            d2Pdy2 = c.*Lx.*d2Ly;
            dPdy   = c.*Lx.*dLy;

            % mixed derivative
            d2Pdxy = c.*dLx.*dLy;

        case 3

            [Lx, dLx] = GradLegendreP(x, i);
            [Ly, dLy] = GradLegendreP(y, j);

            % x-derivative
            dPdx = c.*dLx.*Ly;

            % y-derivative
            dPdy = c.*Lx.*dLy;
    
        case 1

            Lx  = LegendreP(x,i);
            Ly  = LegendreP(y,j);
        
    end

    % psi
    P = c.*Lx.*Ly;

end
