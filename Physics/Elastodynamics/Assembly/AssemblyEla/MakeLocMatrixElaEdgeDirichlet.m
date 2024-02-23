%> @file  MakeLocMatrixElaEdgeDirichlet.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local matrix computation for Dirichlet boundary condition
%>
%==========================================================================
%> @section classMakeLocMatrixElaEdgeDirichlet Class description
%==========================================================================
%> @brief Local matrix computation for Dirichlet boundary condition
%
%> @param Edge        Struct with local edge matrices
%> @param ds          Scaled weight numerical integration
%> @param phiedgeq    Basis function at quadrature nodes
%> @param gradedgeqx  Gradient-x of the basis function at quadrature nodes
%> @param gradedgeqy  Gradient-y of the basis function at quadrature nodes
%> @param par         Physical parameters 
%> @param nx,ny       Normal vector compoenete
%> @param penalty_geom_iedg   Penalty parameters for the current edge
%
%> @retval Edge      Struct with local edge matrices
%>
%==========================================================================

function   [Edge] = MakeLocMatrixElaEdgeDirichlet(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, par, nx, ny, penalty_geom_iedg)


Edge.S1_P_loc = Edge.S1_P_loc + penalty_geom_iedg * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;
Edge.S4_P_loc = Edge.S4_P_loc + penalty_geom_iedg * (ds .* par.harm_ave .* phiedgeq)' * phiedgeq;

Edge.IT1_P_loc = Edge.IT1_P_loc + (ds .* (nx * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqx + ny * par.mu_ave .* gradedgeqy ))' * phiedgeq;
Edge.IT2_P_loc = Edge.IT2_P_loc + (ds .* (ny * par.lambda_ave .* gradedgeqx + nx * par.mu_ave .* gradedgeqy ))' * phiedgeq;
Edge.IT3_P_loc = Edge.IT3_P_loc + (ds .* (ny * par.mu_ave .* gradedgeqx + nx * par.lambda_ave .* gradedgeqy ))' * phiedgeq;
Edge.IT4_P_loc = Edge.IT4_P_loc + (ds .* (nx * par.mu_ave .* gradedgeqx + ny * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqy ))' * phiedgeq;
