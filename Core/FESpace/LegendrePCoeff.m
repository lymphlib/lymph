%> @file  LegendrePCoeff.m
%> @author Mattia Corti
%> @date 6 November 2023 
%> @brief The function constructs the coefficients of the Legendre polynomials
%> 
%> The function constructs the coefficients of the Legendre polynomial of order
%> \f$N\f$. The coefficients of the polynomials are %> hard-coded until the degree 7.
%> Then we compute them by means of the recursive formula \cite handbookformulas :
%> \f[ P_{i+1}(x) = \frac{1}{i+1}\{(2i+1)xP_i(x)-iP_{i-1}(x)\} \f]
%>
%======================================================================
%> @section classDecr Class description
%======================================================================
%> @brief The function evaluates the Legendre polynomials.
%>
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$.
%>
%> @retval P    Matrix containing in the \f$j\f$-th column the coefficients of a Legendre polynomial of order \f$N_j\f$$.
%======================================================================

function [P] = LegendrePCoeff(N)

    % Initialization of the matrix P
    P = zeros(max(N+1), length(N));

    if any(N == 0)
        P(1, N == 0) = 1;
    end

    if any(N == 1)
        P(2, N == 1) = 1;
    end

    if any(N == 2)
        P(3, N == 2) = 3/2;
        P(1, N == 2) = -1/2;
    end

    if any(N == 3)
        P(4, N == 3) = 5/2;
        P(2, N == 3) = -3/2;
    end

    if any(N == 4)
        P(5, N == 4) = 35/8;
        P(3, N == 4) = -30/8;
        P(1, N == 4) = 3/8;
    end

    if any(N == 5)
        P(6, N == 5) = 63/8;
        P(4, N == 5) = -70/8;
        P(2, N == 5) = 15/8;
    end

    if any(N == 6)
        P(7, N == 6) = 231/16;
        P(5, N == 6) = -315/16;
        P(3, N == 6) = 105/16;
        P(1, N == 6) = -5/16;
    end

    if any(N == 7)
        P(8, N == 7) = 429/16;
        P(6, N == 7) = -693/16;
        P(4, N == 7) = 315/16;
        P(2, N == 7) = -35/16;
    end

    P_old  = zeros(max(N+1),1);
    P_oold = zeros(max(N+1),1);

    if any(N > 7)
        % Extraction of the list of polynomials degree larger than 7
        N_app = N;
        N_app(N_app<8) = [];
        N_app = unique(N_app);
     
        % Starting from the last two polynomials
        P_old(8) = 429/16;
        P_old(6) = -693/16;
        P_old(4) = 315/16;
        P_old(2) = -35/16;

        P_oold(7) = 231/16;
        P_oold(5) = -315/16;
        P_oold(3) = 105/16;
        P_oold(1) = -5/16;
    
        start = 7;

        % Cycle over the orders of Legendre polynomials
        for jj = 1:length(N_app)
    
            %Cycle of recursive formula reconstruction
            for ii = start:N_app(jj)-1
                    P_app = 1/(ii+1)*( (2*ii+1)*[0; P_old(1:end-1)] - ii*P_oold);
                    P_oold = P_old;
                    P_old = P_app; 
            end

            % Update the beginning of the next recursive computation
            start = N_app(jj);

            % Update of the matrix
            P(:, N==N_app(jj)) = P_app.*ones(1,sum(N==N_app(jj)));
        end  
    end
end
