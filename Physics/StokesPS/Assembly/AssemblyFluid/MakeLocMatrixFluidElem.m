%> @file  MakeLocMatrixFluidElem.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Local assembly of the matrices for the stokes problem
%>
%==========================================================================
%> @section classMakeLocMatrixFluidElem Class description
%==========================================================================
%> @brief Local assembly of the matrices for the stokes problem
%>
%> @param El   Struct with local matrices
%> @param dx   Differential for integrals
%> @param phiq  Basis functions at quad points
%> @param gradqx  Gradient (x) of basis functions at quad points
%> @param gradqy  Gradient (y) of basis functions at quad points
%> @param par  Struct with physical parameters
%>
%> @retval El Struct with local matrices
%>
%==========================================================================

function  [El] = MakeLocMatrixFluidElem(El, dx, phiq, gradqx, gradqy, par)

%El.MFF_loc   = El.MFF_loc    + (dx .* ( 2*par.mu.^(-1)  .* phiq))' * phiq;
%El.MF_1_loc  = El.MF_1_loc   + (dx .* (-2*par.mu.^(-1)  .* phiq))' * phiq;
El.MF_2_loc  = El.MF_2_loc   + (dx .* (   par.mu.^(-1)  .* phiq))' * phiq;

El.Mprj_loc  = El.Mprj_loc   + (dx .* ( 1               .* phiq))' * phiq;

El.B1_loc    = El.B1_loc + (dx .* gradqx)' * gradqx;
El.B2_loc    = El.B2_loc + (dx .* gradqx)' * gradqy;
El.B3_loc    = El.B3_loc + (dx .* gradqy)' * gradqx;
El.B4_loc    = El.B4_loc + (dx .* gradqy)' * gradqy;

