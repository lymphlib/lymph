%> @file  GauLeg.m
%> @author Mattia Corti
%> @date 3 May 2023
%> @brief Gauss-Legendre nodes and weights on a given interval
%>
%==========================================================================
%> @section classGauLeg Class description
%==========================================================================
%> @brief This routine computes the n Gauss-Legendre nodes and weights
%> on a given interval (a,b)
%>
%> @param a        First interval bound
%> @param b        Second interval bound
%> @param n        Number of Gauss-Legendre nodes to compute 
%>
%> @retval x_GL    Gauss-Legendre nodes
%> @retval w_GL    Gauss-Legendre weights    
%>
%==========================================================================



% This routine computes the n Gauss-Legendre nodes and weights
% on a given interval (a,b)
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [x_GL, w_GL] = GauLeg(a, b, n)

    % Number of the nodes to compute (before reflection)
    m  = (n+1)/2;

    % Computation of the function
    z = cos(pi*((1:m)-0.25)/(n+0.5))';
    
    % Error initialization in scalar and vectorial form
    Err  = eps + 1;
    ErrV = Err*ones(m, 1);

    % Algorithm
    while Err > eps

        p = [ones(m,1) zeros(m,2)]; 
       
        % Iterations
        for j = 1:n
            p(:,2:3) = p(:,1:2);
            p(:,1) = ((2*j-1)/j) * z.*p(:,2) - ((j-1.0)/j) * p(:,3);
        end
 
        % Update of the new step solution
        p_app  = n * (z.*p(:,1)-p(:,2)) ./ (z.^2-1);
        z_old  = z;
        z      = z_old - (p(:,1)./p_app) .* (ErrV > eps);
        
        % Compute the new errors
        ErrV = abs(z-z_old);
        Err = max(ErrV);

    end

    % Computation of the GL nodes in the requested interval
    x_GL  = ((b+a)/2)-((b-a)/2)*[z(1:end-1); -flip(z)];
    
    % Reflection of z and p_app
    z     = [z(1:end-1);     flip(z)];
    p_app = [p_app(1:end-1); flip(p_app)];

    % Computation of the GL weights in the requested interval
    w_GL  = (b-a)./((1-z.^2) .* (p_app.^2));
    
end