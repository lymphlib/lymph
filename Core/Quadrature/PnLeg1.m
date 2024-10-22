%> @file  PnLeg1.m
%> @author Ilario Mazzieri, Paola Gervasio
%> @date 8 October 2024
%> @brief Evaluates the first derivative of Legendre polynomial of degree n
%>
%==========================================================================
%> @section classPnLeg1 Class description
%==========================================================================
%> @brief Evaluates of \f$(L_n)'(x)\f$, \f$L_n(x)\f$ at the node(s) x, using the three 
%> term relation (2.3.3), page 75 \cite CHQZ2.
%>
%> @param x        Scalar or column array of nodes 
%> @param n        Degree of the Legendre polynomial
%>
%> @retval p1    Scalar or column array with the evaluation of (L_n)'(x)
%> @retval p     Scalar or column array with the evaluation of L_n(x)
%>
%==========================================================================

function [p1,p] = PnLeg1 (x, n)
    
    nn = size(x);
    p  = zeros(nn);
    
    if n==0
        p  = ones(nn);
        p1 = zeros(nn);
    elseif n==1
        p  = x;
        p1 = ones(nn);
    else
        p1  = 1 + 0*p;
        p2  = x + 0*p;
        p11 =     0*p;
        p21 = 1 + 0*p;
        p3  =     0*p;
        p31 =     0*p;
        for k = 1 : n-1
            duekp1 = 2*k+1;
            kp1 = 1/(k+1);
            p3  = (duekp1*x.*p2-k*p1)*kp1;
            p31 = (duekp1*(x.*p21+p2)-k*p11)*kp1;
            p11 = p21;
            p21 = p31;
            p1  = p2;
            p2  = p3;
        end
        p1 = p31;
        p  = p3;
    end
end