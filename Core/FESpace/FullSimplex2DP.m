%> @file  FullSimplex2DP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 16 September 2025 
%> @brief Evaluate the 2D modal basis and its derivatives hypercube at 
%> \f$(x,y)\f$ of order \f$(i,j)\f$. The Jacobi polynomial basis is
%> used for triangular meshes, because it is an orthonormal basis on
%> simplices.
%>
%> 
%======================================================================
%> @section classFullSimplex2DP Class description
%======================================================================
%> @brief Evaluate the 2D modal basis and its derivatives hypercube at
%> \f$(x,y)\f$ of order \f$(i,j)\f$. The Jacobi polynomial basis is
%> used for triangular meshes, because it is an orthonormal basis on
%> simplices.
%>
%> @param x     x-coordinates where evaluate the polynomials.
%> @param y     y-coordinates where evaluate the polynomials.
%> @param i     Orders for the construction of the polynomials in x-direction.
%> @param j     Orders for the construction of the polynomials in y-direction.
%>
%> @retval P       2D polynomial on the simplex.
%> @retval dPdr    Derivative in r-direction of 2D polynomial on the simplex.
%> @retval dPds    Derivative in s-direction of 2D polynomial on the simplex.
%======================================================================

function [P, dPdr, dPds] = FullSimplex2DP(x, y, i, j)

    % Computation of 1D Jacobi polynomials
    Jx = JacobiP(x, 0, 0, i);
    Jy = zeros(size(Jx));
    
    if nargout == 3
        % Computation of 1D derivatives of Jacobi polynomials
        dJx = GradJacobiP(x, 0, 0, i);
        dJy = zeros(size(dJx));
    end

    i_unq = unique(i);

    for kk = 1:max(i)+1

        pos = (i==i_unq(kk));
        % Computation of 1D Jacobi polynomials
        Jy(:,pos) = JacobiP(y, 2*i_unq(kk)+1,0, j(pos));

        if nargout == 3
            
            % Computation of 1D derivatives of Jacobi polynomials
            dJy(:,pos) = GradJacobiP(y, 2*i_unq(kk)+1, 0, j(pos));
        end
    end

    if nargout == 3
    
        pos = (i>0);
        % Computation of r-derivative
        dPdr = dJx.*Jy;

        dPdr(:,pos) = dPdr(:,pos).*((0.5*(1-y)).^(i(pos)-1));
    
        % Computation of s-derivative
        dPds = dJx.*(Jy.*(0.5*(1+x)));
    
        dPds(:,pos) = dPds(:,pos).*((0.5*(1-y)).^(i(pos)-1));
        
        tmp = dJy.*((0.5*(1-y)).^i);

        tmp(:,pos) = tmp(:,pos) - 0.5*i(pos).*Jy(:,pos).*((0.5*(1-y)).^(i(pos)-1));
        
        dPds = dPds + Jx.*tmp;

        % Normalization of the final gradients
        dPdr = 2*2.^(i+0.5).*dPdr;
        dPds = 2*2.^(i+0.5).*dPds;

    end

    % Computation of 2D Jacobi polynomials
    P = sqrt(2.0)*Jx.*Jy.*((1-y).^i);

end
