%> @file  GauLeg.m
%> @author Ilario Mazzieri, Paola Gervasio, Mattia Corti
%> @date 8 October 2024
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
%> @param n        Number of Gauss-Legendre nodes to compute (n odd, n>1)
%>
%> @retval x_GL    Gauss-Legendre nodes
%> @retval w_GL    Gauss-Legendre weights
%>
%==========================================================================

function [x_GL, w_GL] = GauLeg(a, b, n)

    if n <= 1
        x_GL = 0;
        w_GL = 2;
    else 
        x = JacobiRoots(n,0,0);
        w = 2./(PnLeg1(x,n).^2.*(1-x.^2));
        
        % Computation of the GL nodes in the requested interval
        x_GL = 0.5*(b-a)*x + 0.5*(b+a);
        w_GL = w * 0.5*(b-a);
    end

end


