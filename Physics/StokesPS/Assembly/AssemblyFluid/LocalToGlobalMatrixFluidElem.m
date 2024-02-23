%> @file  LocalToGlobalMatrixFluidElem.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Global assembly of the matrices for the stokes problem
%>
%==========================================================================
%> @section classLocalToGlobalMatrixFluidElem Class description
%==========================================================================
%> @brief Global assembly of the matrices for the stokes problem
%>
%> @param A   Struct with Global matrices
%> @param El Struct with local matrices
%> @param index  mapping from local to global indeces
%>
%> @retval A   Struct with Global matrices
%>
%==========================================================================

function [A] = LocalToGlobalMatrixFluidElem(A, El, index)

A.MFF(index,index)   = El.MFF_loc;
A.MF_1(index,index)  = El.MF_1_loc;
A.MF_2(index,index)  = El.MF_2_loc;
A.Mprj(index,index)  = El.Mprj_loc;

A.B1(index,index)    = El.B1_loc;
A.B2(index,index)    = El.B2_loc;
A.B3(index,index)    = El.B3_loc;
A.B4(index,index)    = El.B4_loc;