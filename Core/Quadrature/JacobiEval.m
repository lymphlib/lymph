%> @file  JacobiEval.m
%> @author Ilario Mazzieri, Paola Gervasio
%> @date 8 October 2024
%> @brief Evaluates Jacobi polynomials
%>
%==========================================================================
%> @section classJacobiEval Class description
%==========================================================================
%> @brief Evaluates Jacobi polynomial \f$P_n^{(\alpha,\beta)}\f$ at \f$x
%> \in(-1,1)\f$ (formula (2.5.3) page 92, \cite CHQZ2).
%>
%> @param x        Scalar or one-dimensional array of length (m)
%> @param n        Polynomial degree
%> @param alpha    First parameter of the Jacobi polynomial
%> @param beta     Second parameter of the Jacobi polynomial
%>
%> @retval p       p = [P_n^{(\alpha,\beta)}(x),
%                       P_(n-1)^{(\alpha,\beta)}(x),
%                       P_(n-2)^{(\alpha,\beta)}(x)];
%
%  @retval pd      pd = [(P_n^{(\alpha,\beta)})'(x),
%                        (P_(n-1)^{(\alpha,\beta)})'(x),
%                        (P_(n-2)^{(\alpha,\beta)})'(x)];
%>
%==========================================================================

function [p,pd] = JacobiEval(x,n,alpha,beta)

    apb = alpha + beta;
    ab2 = alpha^2 - beta^2;
    
    if isscalar(x)
        % x is a scalar
        p  = 1;  pd  = 0;
        p1 = 0;  pd1 = 0;
        p2 = 0;  pd2 = 0;
        if n == 0
            return
        elseif n == 1
            p1 = p; 
            p2 = p1;
            p = (alpha-beta+(apb+2)*x)*0.5;
    
            pd1 = pd; 
            pd2 = pd1;
            pd = 0.5*(apb+2);
        else
            p1 = p; 
            p2 = p1;
            p = (alpha-beta+(apb+2)*x)*0.5;
    
            pd1 = pd; 
            pd2 = pd1;
            pd = 0.5*(apb+2);
    
            for k = 1:n-1
                k1 = k+1; 
                k2 = k*2; 
                k2ab  = k2 + alpha + beta;
                k2ab1 = k2ab + 1; 
                k2ab2 = k2ab1 + 1;
                p2 = p1; 
                p1 = p;
                pd2 = pd1; 
                pd1 = pd;
                a1 = 2*k1*(k1+apb)*k2ab;
                % Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x
                a21 = k2ab1*ab2;
                a22 = k2ab2*k2ab1*k2ab;
                a3 = 2*(k+alpha)*(k+beta)*k2ab2;
                p  = ((a21+a22*x)*p1-a3*p2)/a1;
                pd = (a22*p1+(a21+a22*x)*pd1-a3*pd2)/a1;
            end
        end
    
    else
        [m1,m2] = size(x);
        if m1 < m2
            x = x';
        end
        m = max(m1,m2);
        p = [ones(m,1),zeros(m,1),zeros(m,1)];   pd=zeros(m,3);
        if n == 0
            return
    
        elseif n==1
            p(:,2) = p(:,1);
            p(:,3) = p(:,2);
            p(:,1) = (alpha-beta+(apb+2)*x)*0.5;
    
            pd(:,2) = pd(:,1);
            pd(:,3) = pd(:,2);
            pd(:,1) = 0.5*(apb+2);
    
        else
            p(:,2) = p(:,1);
            p(:,3) = p(:,2);
            p(:,1) = (alpha-beta+(apb+2)*x)*0.5;
    
            pd(:,2) = pd(:,1);
            pd(:,3) = pd(:,2);
            pd(:,1) = 0.5*(apb+2);
    
            for k = 1 : n-1
                k2 = k*2;
                k2ab = k2 + alpha + beta;
                p(:,3)  = p(:,2);
                p(:,2)  = p(:,1);
                pd(:,3) = pd(:,2);
                pd(:,2) = pd(:,1);
                a1 = 2*(k+1)*(k+apb+1)*k2ab;
                % Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x
                a21 = (k2ab+1)*ab2;
                a22 = (k2ab+2)*(k2ab+1)*k2ab;
                a3  = 2*(k+alpha)*(k+beta)*(k2ab+2);
                p(:,1) = ((a21+a22*x).*p(:,2)-a3*p(:,3))/a1;
                pd(:,1)= (a22*p(:,2)+(a21+a22*x).*pd(:,2)-a3*pd(:,3))/a1;
            end
        end
    
    end

end
