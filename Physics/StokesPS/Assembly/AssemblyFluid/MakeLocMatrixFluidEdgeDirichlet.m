%> @file  MakeLocMatrixFluidEdgeDirichlet.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local matrix computation for Dirichlet boundary condition
%>
%==========================================================================
%> @section classMakeLocMatrixFluidEdgeDirichlet Class description
%==========================================================================
%> @brief Local matrix computation for Dirichlet boundary condition
%>
%> @param Edge        Struct with local edge matrices
%> @param ds          Scaled weight numerical integration
%> @param phiedgeq    Basis function at quadrature nodes
%> @param gradedgeqx  Gradient-x of the basis function at quadrature nodes
%> @param gradedgeqy  Gradient-y of the basis function at quadrature nodes
%> @param nx,ny       Normal vector compoenete
%> @param penalty_geom_iedg   Penalty parameters for the current edge
%> 
%> @retval Edge      Struct with local edge matrices
%>
%==========================================================================

function   [Edge] = MakeLocMatrixFluidEdgeDirichlet(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, nx, ny, penalty_geom_iedg)


Edge.S1_B_loc = Edge.S1_B_loc + penalty_geom_iedg * (ds .* phiedgeq)' * phiedgeq * nx * nx;
Edge.S2_B_loc = Edge.S2_B_loc + penalty_geom_iedg * (ds .* phiedgeq)' * phiedgeq * nx * ny;
Edge.S3_B_loc = Edge.S3_B_loc + penalty_geom_iedg * (ds .* phiedgeq)' * phiedgeq * ny * nx;
Edge.S4_B_loc = Edge.S4_B_loc + penalty_geom_iedg * (ds .* phiedgeq)' * phiedgeq * ny * ny;

Edge.BT1_loc = Edge.BT1_loc + (ds .* ( nx * gradedgeqx))' * phiedgeq;
Edge.BT2_loc = Edge.BT2_loc + (ds .* ( ny * gradedgeqx))' * phiedgeq;
Edge.BT3_loc = Edge.BT3_loc + (ds .* ( nx * gradedgeqy))' * phiedgeq;
Edge.BT4_loc = Edge.BT4_loc + (ds .* ( ny * gradedgeqy))' * phiedgeq;
