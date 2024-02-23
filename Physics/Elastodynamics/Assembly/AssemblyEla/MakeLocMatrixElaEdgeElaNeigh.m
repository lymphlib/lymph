%> @file  MakeLocMatrixElaEdgeElaNeigh.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local matrix computation for elastic element (edge matrix for
%neighboring elements)
%>
%==========================================================================
%> @section classMakeLocMatrixElaEdgeElaNeigh Class description
%==========================================================================
%> @brief Local matrix computation for elastic element (edge matrix for
%neighboring elements)
%
%> @param Edge        Struct with local edge matrices
%> @param ds          Scaled weight numerical integration
%> @param phiedgeq    Basis function at quadrature nodes
%> @param gradedgeqx  Gradient-x of the basis function at quadrature nodes
%> @param gradedgeqy  Gradient-y of the basis function at quadrature nodes
%> @param phiedgeqneigh    Basis function of neigh el at quadrature nodes
%> @param par         Physical parameters 
%> @param nx,ny       Normal vector compoenete
%> @param penalty_geom_iedg   Penalty parameters for the current edge
%> @param iedg        Index of current edge
%
%> @retval Edge      Struct with local edge matrices
%>
%==========================================================================

function  [Edge] = MakeLocMatrixElaEdgeElaNeigh(Edge, ds, phiedgeq, gradedgeqx, gradedgeqy, phiedgeqneigh, par, nx, ny, penalty_geom_iedg, iedg)

Edge.SN1_loc(:,:,iedg) = Edge.SN1_loc(:,:,iedg) - (ds .* (par.harm_ave .* penalty_geom_iedg) .* phiedgeq)' * phiedgeqneigh;
Edge.SN4_loc(:,:,iedg) = Edge.SN4_loc(:,:,iedg) - (ds .* (par.harm_ave .* penalty_geom_iedg) .* phiedgeq)' * phiedgeqneigh;

Edge.ITN1_loc(:,:,iedg) = Edge.ITN1_loc(:,:,iedg) - 0.5 * (ds .* (nx * (par.lambda_ave + 2*par.mu_ave).* gradedgeqx + ny * par.mu_ave .* gradedgeqy ))' * phiedgeqneigh;
Edge.ITN2_loc(:,:,iedg) = Edge.ITN2_loc(:,:,iedg) - 0.5 * (ds .* (ny * par.lambda_ave .* gradedgeqx + nx * par.mu_ave .* gradedgeqy ))' * phiedgeqneigh;
Edge.ITN3_loc(:,:,iedg) = Edge.ITN3_loc(:,:,iedg) - 0.5 * (ds .* (ny * par.mu_ave .* gradedgeqx + nx * par.lambda_ave .* gradedgeqy ))' * phiedgeqneigh;
Edge.ITN4_loc(:,:,iedg) = Edge.ITN4_loc(:,:,iedg) - 0.5 * (ds .* (nx * par.mu_ave .* gradedgeqx + ny * (par.lambda_ave + 2*par.mu_ave) .* gradedgeqy ))' * phiedgeqneigh;

