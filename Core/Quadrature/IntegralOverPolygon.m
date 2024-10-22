%> @file  IntegralOverPolygon.m
%> @author Giorgio Pennesi, Mattia Corti, Paul Houston
%> @date 11 November 2023 
%> @brief Evaluate integral of a polygon of verteces v
%> 
%> The function computes the integrals of all the two-dimensional monomials
%> until the degree ex in the x-component and ey in the y-component and
%> with coefficients Cf over the polygon of vertices v. Namely:
%> \f[
%>      f(x,y) = Cf x^i y^j,\qquad i=1,...,E_x,\;j=1,...,E_y
%> \f]
%> For the details of formulas, we refer to \cite pennesiquadraturefree
%>
%======================================================================
%> @section classIntegralOverPolygon Class description
%======================================================================
%> @brief Evaluate integral of a polygon of verteces v.
%>
%> @param Cf        Multiplicative coefficients of the monomials.
%> @param ex        Maximum degree of the monomials in x-component
%> @param ey        Maximum degree of the monomials in y-component
%> @param v         Vertices of the polygon
%>
%> @retval I    Computed volume integral associated to the 2D-monomials.
%> @retval If   Computed surface integrals associated to the 2D-monomials.
%======================================================================

function [I, If] = IntegralOverPolygon(Cf, ex, ey, v)

    % Space dimension
    d = size(v,1);
    
    % Number of vertices of the polygon
    m = size(v,2);
    
    % Initialization of the structures
    I = zeros(ex+1,ey+1);
    If = zeros(ex+1,ey+1,m);
    Ex = (0:ex)';
    Ey = 0:ey;
    
    % Cycle over the polygon vertices
    for ii = 1:m
        
        % Initialization of the structures
        FHM = zeros(ex+1,ey+1);
        var_x = ones(1,ex);
        var_y = ones(1,ey);

        % Extraction of the first point coordinates
        x1 = v(1,ii);        y1 = v(2,ii);

        % Extraction of the second point coordinates
        if ii < m
            x2 = v(1,ii+1);   y2 = v(2,ii+1);
        else
            x2 = v(1,1);      y2 = v(2,1);
        end
        
        % Construction of the parameter bi and dij
        bi = ((y2-y1)*x1 + (x1-x2)*y1) / sqrt( (y2-y1)^2 + (x1-x2)^2 );
        dij = norm([x2-x1,y2-y1]);
        
        % Integral of the constant on the edge
        FHM(1,1) = dij/(d-1);
        
        % Integral of the monomials only in x

        for k=1:ex
            if k==1
                var_x(1) = x2;
            else
                var_x(k) = var_x(k-1)*x2;
            end
            FHM(k+1,1) = ( dij*var_x(k) + x1*k*FHM(k,1) )/(d+k-1);
        end
        

        % Integral of the monomials only in y
        for k=1:ey
            if k==1
                var_y(1) = y2;
            else
                var_y(k) = var_y(k-1)*y2;
            end
            FHM(1,k+1) = ( dij*var_y(k) + y1*k*FHM(1,k) )/(d+k-1);
        end
        
        % Integral of the monomials in x and y
        for k=1:ex
            for l=1:ey
                FHM(k+1,l+1) = ( dij*var_x(k)*var_y(l) + x1*k*FHM(k,l+1) + y1*l*FHM(k+1,l) )/(d+k+l-1);
            end
        end
        
        % Final integrals computation of the polynomials basis
        I = I + bi*FHM./(d+Ex+Ey);
        % Save the faces integrals in a dedicated structure
        If(:,:,ii) = FHM;
    end

    % Final integral computation
    I = Cf*I;

    If = Cf*If;

    % Reshape for optimization purposes
    I = reshape(I,[size(I,1)*size(I,2) 1]);
    If = reshape(If,[size(If,1)*size(If,2) m]);
end
