%> @file  GradSimplex2DP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief Evaluate 2D derivatives of the modal basis on simplex at \f$(a,b)\f$ of order \f$(i,j)\f$
%> 
%> The function evaluates 2D derivatives of the modal basis on simplex at 
%> \f$(a,b)\f$ of order \f$(i,j)\f$.
%>
%> 
%======================================================================
%> @section classGradSimplex2DP Class description
%======================================================================
%> @brief Evaluate 2D derivatives of the modal basis on simplex at
%> \f$(a,b)\f$ of order \f$(i,j)\f$.
%>
%>
%> @param a     x-coordinates where evaluate the polynomials.
%> @param b     y-coordinates where evaluate the polynomials.
%> @param i     Orders for the construction of the polynomials in x-direction.
%> @param j     Orders for the construction of the polynomials in y-direction.
%>
%> @retval dPdr    Derivative in r-direction of 2D polynomial on the simplex.
%> @retval dPds    Derivative in s-direction of 2D polynomial on the simplex.
%======================================================================

function [dPdr, dPds] = GradSimplex2DP(a, b, i, j)

    % Computation of 1D Jacobi polynomials
    Fx = JacobiP(a, 0, 0, i);
    Gy = JacobiP(b, 2*i+1,0, j);

    % Computation of 1D derivatives of Jacobi polynomials
    dFx = GradJacobiP(a, 0, 0, i);
    dGy = GradJacobiP(b, 2*i+1,0, j);

    % Computation of r-derivative: d/dr = da/dr d/da + db/dr d/db = (2/(1-b)) d/da
    dPdr = dFx.*Gy;
    
    if(i>0)
        dPdr = dPdr.*((0.5*(1-b)).^(i-1));
    end
    
    % Computation of s-derivative: d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
    dPds = dFx.*(Gy.*(0.5*(1+a)));
    
    if(i>0)
        dPds = dPds.*((0.5*(1-b)).^(i-1));
    end

    tmp = dGy.*((0.5*(1-b)).^i);

    if(i>0)
        tmp = tmp - 0.5*i*Gy.*((0.5*(1-b)).^(i-1));
    end

    dPds = dPds + Fx.*tmp;

    % Normalization of the final gradients
    dPdr = 2*2^(i+0.5)*dPdr;
    dPds = 2*2^(i+0.5)*dPds;

end
