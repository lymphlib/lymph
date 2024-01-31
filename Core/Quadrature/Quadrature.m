%> @file  Quadrature.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 3 May 2023
%> @brief Computation of the n Gauss-Legendre nodes and weights.
%>
%==========================================================================
%> @section classQuadrature Class description
%==========================================================================
%> @brief This routine computes the n Gauss-Legendre nodes and weights on the
%>  reference interval (0,1)  to be used for face integrals, and
%>  the n^2 Gauss-Legendre nodes and weights on the reference triangle
%>  (0,0), (1,0), (0,1) or on the reference square (-1,1)x(-1,1) to be used
%>  for volume integrals.
%>
%> @param n         Number of Gauss-Legendre nodes to compute
%>
%> @retval node_1D     Gauss-Legendre nodes on the reference interval (0,1)
%> @retval w_1D        Gauss-Legendre weights on the reference interval (0,1)
%> @retval node_2D     Gauss-Legendre nodes on the reference interval (-1,1)
%> @retval w_2D        Gauss-Legendre weights on the reference interval (-1,1)
%>
%==========================================================================

function [node_1D, w_1D, node_2D, w_2D] = Quadrature(n)

    % 1D Gauss-Legendre nodes on the interval [0 1]
    [node_1D, w_1D] = GauLeg(0,1,n);

    % 1D Gauss-Legendre nodes on the interval [-1 1]
    [x_GL, w_GL] = GauLeg(-1,1,n);

    X1 = reshape((ones(n,n).*x_GL'), [n^2 1]);
    X2 = reshape((ones(n,n).*x_GL),  [n^2 1]);
    W1 = reshape((ones(n,n).*w_GL'), [n^2 1]);
    W2 = reshape((ones(n,n).*w_GL),  [n^2 1]);

    % 2D Gauss-Legendre nodes on the set [-1 1]x[-1 1]
    node_2D = [(1+X1)/2 (1-X1).*(1+X2)/4];
    w_2D    = (1-X1).*W1.*W2/8;

end

