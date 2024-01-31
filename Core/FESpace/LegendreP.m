%> @file  LegendreP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief The function evaluates the Legendre polynomials
%> 
%> The function evaluates the Legendre polynomial of order
%> \f$N\f$ at point \f$x\f$, namely \f$P_N(x)\f$. The polynomials are
%> hard-coded until the degree 7. Then we compute them by means of the recursive formula \cite handbookformulas :
%> \f[ P_{i+1}(x) = \frac{1}{i+1}\{(2i+1)xP_i(x)-iP_{i-1}(x)\} \f]
%>
%======================================================================
%> @section classLegendreP Class description
%======================================================================
%> @brief The function evaluates the Legendre polynomials.
%>
%> @param x     Vector of points \f$x_i\f$ in which the Legendre polynomials are evaluated.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval P    Matrix containing at the position \f$P_{ij}\f$ a Legendre polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$.
%======================================================================

function [P] = LegendreP(x, N)
    
    % Turn points into column if needed
    x = x(:); 

    % Initialization of the matrix P
    P = zeros(length(x), length(N));

    if any(N == 0)
        P = P + (N == 0).*ones(size(x));
    end

    if any(N == 1)
        P = P + (N == 1).*x;
    end

    if any(N == 2)
        P = P + (N == 2).*(3*x.^2-1)/2;
    end

    if any(N == 3)
        P = P + (N == 3).*(5*x.^3-3*x)/2;
    end

    if any(N == 4)
        P = P + (N == 4).*(35*x.^4-30*x.^2+3)/8;
    end

    if any(N == 5)
        P = P + (N == 5).*(63*x.^5-70*x.^3+15*x)/8;
    end

    if any(N == 6)
        P = P + (N == 6).*(231*x.^6-315*x.^4+105*x.^2-5)/16;
    end

    if any(N == 7)
        P = P + (N == 7).*(429*x.^7-693*x.^5+315*x.^3-35.*x)/16;
    end

    if any(N > 7)
        % Extraction of the list of polynomials degree larger than 7
        N_app = N;
        N_app(N_app<8) = [];
        N_app = unique(N_app);
        
        % Starting from the last two polynomials
        P_old  = (429*x.^7-693*x.^5+315*x.^3-35.*x)/16;
        P_oold = (231*x.^6-315*x.^4+105*x.^2-5)/16;
        start = 7;

        % Cycle over the orders of Legendre polynomials
        for jj = 1:length(N_app)
            
            % Cycle of recursive formula reconstruction
            for ii = start:N_app(jj)-1
                    P_app = 1/(ii+1)*( (2*ii+1)*x.*P_old - ii*P_oold);
                    P_oold = P_old;
                    P_old = P_app; 
            end
    
            % Update the beginning of the next recursive computation
            start = N_app(jj);

            % Update of the matrix
            P = P + (N==N_app(jj)).*P_app;
        end  
    end
end
