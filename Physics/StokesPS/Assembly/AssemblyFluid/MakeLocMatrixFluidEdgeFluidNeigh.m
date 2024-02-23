%> @file  MakeLocMatrixFluidEdgeFluidNeigh.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local matrix computation for fluid element (edge matrix for
%neighboring elements)
%>
%==========================================================================
%> @section classMakeLocMatrixFluidEdgeFluidNeigh Class description
%==========================================================================
%> @brief Local matrix computation for fluid element (edge matrix for
%neighboring elements)
%>
%> @param Edge        Struct with local edge matrices
%> @param ds          Scaled weight numerical integration
%> @param phiedgeq    Basis function at quadrature nodes
%> @param gradedgeqx  Gradient-x of the basis function at quadrature nodes
%> @param gradedgeqy  Gradient-y of the basis function at quadrature nodes
%> @param phiedgeqneigh    Basis function of neigh el at quadrature nodes
%> @param nx,ny       Normal vector compoenete
%> @param penalty_geom_iedg   Penalty parameters for the current edge
%> @param iedg        Index of current edge
%> 
%> @retval Edge      Struct with local edge matrices
%>
%==========================================================================

function  [Edge] = MakeLocMatrixFluidEdgeFluidNeigh(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, phiedgeqneigh, nx, ny, penalty_geom_iedg, iedg)

Edge.SN_B1_loc(:,:,iedg) = Edge.SN_B1_loc(:,:,iedg) - (ds .* ( penalty_geom_iedg) .* phiedgeq)' * phiedgeqneigh * nx * nx;
Edge.SN_B2_loc(:,:,iedg) = Edge.SN_B2_loc(:,:,iedg) - (ds .* ( penalty_geom_iedg) .* phiedgeq)' * phiedgeqneigh * nx * ny;
Edge.SN_B3_loc(:,:,iedg) = Edge.SN_B3_loc(:,:,iedg) - (ds .* ( penalty_geom_iedg) .* phiedgeq)' * phiedgeqneigh * ny * nx;
Edge.SN_B4_loc(:,:,iedg) = Edge.SN_B4_loc(:,:,iedg) - (ds .* ( penalty_geom_iedg) .* phiedgeq)' * phiedgeqneigh * ny * ny;

Edge.BTN1_loc(:,:,iedg) = Edge.BTN1_loc(:,:,iedg) - 0.5 * (ds .*  ( nx * gradedgeqx))' * phiedgeqneigh;
Edge.BTN2_loc(:,:,iedg) = Edge.BTN2_loc(:,:,iedg) - 0.5 * (ds .*  ( ny * gradedgeqx))' * phiedgeqneigh;
Edge.BTN3_loc(:,:,iedg) = Edge.BTN3_loc(:,:,iedg) - 0.5 * (ds .*  ( nx * gradedgeqy))' * phiedgeqneigh;
Edge.BTN4_loc(:,:,iedg) = Edge.BTN4_loc(:,:,iedg) - 0.5 * (ds .*  ( ny * gradedgeqy))' * phiedgeqneigh;