%> @file  JacobiP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief The function evaluates the Jacobi polynomials
%> 
%> The function evaluates the Jacobi polynomial with parameters \f$\alpha\f$ and \f$\beta\f$ 
%> of order \f$N\f$ at point \f$x\f$, namely \f$P_N^{(\alpha,\beta)}(x)\f$. The polynomials are
%> hard-coded until the degree 1. Then we compute them by means of the recursive formula \cite jacobiformulas.
%>

%> Note: They are normalized to be orthonormal.
%> 
%======================================================================
%> @section classJacobiP Class description
%======================================================================
%> @brief The function evaluates the Jacobi polynomials.
%>
%> @param x     Vector of points \f$x_i\f$ in which the Jacobi polynomials are evaluated.
%> @param alpha First numerical parameter of Jacobi polynomials.
%> @param beta  Second numerical parameter of Jacobi polynomials.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval P   Matrix containing at the position \f$P_{ij}\f$ a Jacobi polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$, computed with parameters \f$(\alpha,\beta)\f$.
%===================================================s===================

function [P] = JacobiP(x, alpha, beta, N)

    % Turn points into row if needed
    x = x(:); 

    % Initialization of the matrix P
    P = zeros(length(x), length(N));

    % Polynomial constants
    gamma0 = (2^(alpha+beta+1)/(alpha+beta+1)) * ... 
        gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+1);

    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;

    if any(N == 0)
        P = P + (N == 0).*ones(size(x))/sqrt(gamma0);
    end

    if any(N == 1)
        P = P + (N == 1).*(((alpha+beta+2)*x/2 + (alpha-beta)/2)/sqrt(gamma1));
    end

    if any(N > 1)
        % Extraction of the list of polynomials degree larger than 7
        N_app = N;
        N_app(N_app<2) = [];
        N_app = unique(N_app);
        
        % Starting from the last two polynomials
        P_old  = (((alpha+beta+2)*x/2 + (alpha-beta)/2)/sqrt(gamma1));
        P_oold = ones(size(x))/sqrt(gamma0);
        start = 1;

        % Cycle over the orders of Legendre polynomials
        for jj = 1:length(N_app)
            
            c_oold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

            % Cycle of recursive formula reconstruction
            for ii = start:N_app(jj)-1
                    
                    % Compute the constants of the recursive formula
                    c_app = 2*ii+alpha+beta;
                    c_new = 2/(c_app+2)*sqrt((ii+1)*(ii+1+alpha+beta)* ...
                        (ii+1+alpha)*(ii+1+beta)/(c_app+1)/(c_app+3));
                    c_old = - (alpha^2-beta^2)/(c_app*(c_app+2));

                    % Recursive step 
                    P_app = 1/c_new*((x-c_old).*P_old - c_oold.*P_oold);
                    
                    % Update of the values
                    P_oold = P_old;
                    P_old = P_app;
                    c_oold = c_new;
            end
    
            % Update the beginning of the next recursive computation
            start = N_app(jj);

            % Update of the matrix
            P = P + (N==N_app(jj)).*P_app;
        end  
    end
end