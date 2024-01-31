%> @file  MakeLocMatrixElaElem.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local assembly of the matrices for the elastic problem \cite AM2017
%>
%==========================================================================
%> @section classMakeLocMatrixElaElem Class description
%==========================================================================
%> @brief Local assembly of the matrices for the elastic problem \cite AM2017
%
%> @param El   Struct with local matrices
%> @param dx   Differential for integrals
%> @param phiq  Basis functions at quad points
%> @param gradqx  Gradient (x) of basis functions at quad points
%> @param gradqy  Gradient (y) of basis functions at quad points
%> @param par  Struct with physical parameters
%
%> @retval El Struct with local matrices 
%>
%==========================================================================

function [El] = MakeLocMatrixElaElem(El, dx, phiq, gradqx, gradqy, par)

% inner index is i - outer index is j
El.V1_loc        = El.V1_loc + (dx .* ((par.lam+2*par.mu) .* gradqx))' * gradqx + (dx .* (par.mu .* gradqy))' * gradqy;
El.V2_loc        = El.V2_loc + (dx .* (par.lam        .* gradqx))' * gradqy + (dx .* (par.mu .* gradqy))' * gradqx;
El.V3_loc        = El.V3_loc + (dx .* (par.lam        .* gradqy))' * gradqx + (dx .* (par.mu .* gradqx))' * gradqy;
El.V4_loc        = El.V4_loc + (dx .* ((par.lam+2*par.mu) .* gradqy))' * gradqy + (dx .* (par.mu .* gradqx))' * gradqx;


El.M1_P_rho_loc  = El.M1_P_rho_loc  + (dx .* (par.rho_e       .* phiq))' * phiq;
El.MPrjP_1_loc   = El.MPrjP_1_loc   + (dx .* (1         .* phiq))' * phiq;
                  
El.D1_loc  = El.D1_loc  + (dx .* (2 * par.rho_e .* par.zeta .* phiq))' * phiq;
El.C1_loc  = El.C1_loc  + (dx .* (par.rho_e .* par.zeta.^2  .* phiq))' * phiq;

