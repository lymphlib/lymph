%> @file  JacobiRoots.m
%> @author Ilario Mazzieri, Paola Gervasio
%> @date 8 October 2024
%> @brief Computes the n zeros of the Jacoby polynomial.
%>
%==========================================================================
%> @section classJacobiRoots Class description
%==========================================================================
%> @brief Computes the n zeros of the Jacoby polynomial 
%> \f$P_n^{(\alpha,\beta)}(x)\f$ by Newton method and deflation process.
%>
%> @param n        Polynomial degree
%> @param alpha    First parameter of the Jacobi polynomial
%> @param beta     Second parameter of the Jacobi polynomial
%>
%> @retval x       Zeros of Jacobi polynomial
%>
%==========================================================================

function [x] = JacobiRoots(n,alpha,beta)

    x    = zeros(n,1);
    x0   = cos(pi/(2*n));
    tol  = 1.e-14;
    kmax = 15;
    
    for j = 1 : n
        diff = tol + 1;
        kiter = 0;
        while kiter <= kmax && diff >= tol
            [p,pd] = JacobiEval(x0,n,alpha,beta);
            % deflation process
            % q(x)       = p(x)*(x-x_1)*... (x-x_{j-1})
            % q(x)/q'(x) = p(x)/[p'(x)-p(x)*\sum_{i<j} 1/(x-x_i)]
            ss    = sum(1./(x0-x(1:j-1)));
            x1    = x0-p/(pd-ss*p);
            diff  = abs(x1-x0);
            kiter = kiter+1;
            x0 = x1;
        end
        x0 = 5.d-1*(x1+cos((2*(j+1)-1)*pi/(2*n)));
        x(j) = x1;
    end
    x = sort(x);

end
