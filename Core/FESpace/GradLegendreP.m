%> @file  GradLegendreP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief The function evaluates the gradient of Legendre polynomial 
%> 
%> The function evaluates the gradient of Legendre polynomial of order
%> \f$N\f$ at point \f$x\f$, namely \f$P_N'(x)\f$. The polynomials are
%> hard-coded until the degree 7. Then we compute them by means of the recursive formula \cite handbookformulas :
%> \f[ P'_{i+1}(x) = P'_{i-1}(x) + (2i+1)P_{i}(x) \f]
%>
%======================================================================
%> @section classGradLegendreP Class description
%======================================================================
%> @brief The function evaluates the gradient of Legendre polynomial.
%>
%>
%> @param x     Vector of points \f$x_i\f$ in which the Legendre polynomials are evaluated.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval dP   Matrix containing at the position \f$P_{ij}'\f$ a gradient of a Legendre polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$.
%======================================================================

function [dP] = GradLegendreP(x,N)
    
    % Turn points into column if needed
    x = x(:); 

    % Initialization of the matrix P equivalent to dP_0
    dP = zeros(length(x), length(N));

    if any(N == 1)
            dP = dP + (N == 1).*ones(size(x));
    end

    if any(N == 2)
            dP = dP + (N==2).*(3*x);
    end    

    if any(N == 3)
            dP = dP + (N==3).*(15*x.^2-3*ones(size(x)))/2;
    end

    if any(N == 4)
            dP = dP + (N==4).*(35*x.^3-15*x)/2;
    end   

    if any(N == 5)
            dP = dP + (N==5).*(315*x.^4-210*x.^2+15*ones(size(x)))/8;
    end    

    if any(N == 6)
            dP = dP + (N==6).*(1386*x.^5-1260*x.^3+210*x)/16;
    end

    if any(N == 7)
            dP = dP + (N==7).*(3003*x.^6-3465*x.^4+945*x.^2-35*ones(size(x)))/16;
    end 

    % Recoursive computation of the remaining polynomials
    if any(N > 7)
        % Extraction of the list of polynomials degree larger than 7
        N_app = N;
        N_app(N_app<8) = [];
        N_app = unique(N_app);
        
        % Cycle over the orders of Legendre polynomials
        for jj = 1:length(N_app)
            
            % Vectors of degrees of polynomials
            ii = N_app(jj)-1:-2:0;
            
            % Computation of necessary Legendre polynomials
            [P] = LegendreP(x,ii);

            % Computation of the derivative of Legendre polynomials
            dP_app = sum((2*ii+1).*P,2);
            dP = dP + (N==N_app(jj)).*dP_app;

        end
    end
end
