%> @file  MakeLocMatrixElaEdgeAbsorbing.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local matrix computation for absorbing conditions
%>
%==========================================================================
%> @section classMakeLocMatrixElaEdgeAbsorbing Class description
%==========================================================================
%> @brief Local matrix computation for absorbing conditions
%
%> @param Edge      Struct with local edge matrices
%> @param ds        Scaled weight numerical integration
%> @param phiedgeq  Basis function at quadrature nodes
%> @param gradedgeqx  Gradient-x of Basis function at quadrature nodes
%> @param gradedgeqy  Gradient-y of Basis function at quadrature nodes
%> @param par       physical parameters 
%
%> @param Edge      Struct with local edge matrices
%>
%==========================================================================

function  [Edge] = MakeLocMatrixElaEdgeAbsorbing(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, par)

Edge.ABC_S1_loc = Edge.ABC_S1_loc + (ds .* par.c11 .* phiedgeq)' * phiedgeq;
Edge.ABC_S2_loc = Edge.ABC_S2_loc + (ds .* par.c12 .* phiedgeq)' * phiedgeq;
Edge.ABC_S3_loc = Edge.ABC_S3_loc + (ds .* par.c21 .* phiedgeq)' * phiedgeq;
Edge.ABC_S4_loc = Edge.ABC_S4_loc + (ds .* par.c22 .* phiedgeq)' * phiedgeq;

Edge.ABC_R1_loc = Edge.ABC_R1_loc + (ds .* (par.c3_11 .* gradedgeqx + par.c4_11 .* gradedgeqy))' * phiedgeq;
Edge.ABC_R2_loc = Edge.ABC_R2_loc + (ds .* (par.c3_12 .* gradedgeqx + par.c4_12 .* gradedgeqy))' * phiedgeq;
Edge.ABC_R3_loc = Edge.ABC_R3_loc + (ds .* (par.c3_21 .* gradedgeqx + par.c4_21 .* gradedgeqy))' * phiedgeq;
Edge.ABC_R4_loc = Edge.ABC_R4_loc + (ds .* (par.c3_22 .* gradedgeqx + par.c4_22 .* gradedgeqy))' * phiedgeq;
