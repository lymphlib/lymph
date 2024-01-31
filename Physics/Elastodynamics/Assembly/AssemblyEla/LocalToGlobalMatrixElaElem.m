%> @file  LocalToGlobalMatrixElaElem.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Global assembly of the matrices for the elastic problem \cite AM2017
%>
%==========================================================================
%> @section classLocalToGlobalMatrixElaElem Class description
%==========================================================================
%> @brief Global assembly of the matrices for the elastic problem \cite AM2017
%
%> @param A   Struct with Global matrices
%> @retval El Struct with local matrices 
%> @param index  mapping from local to global indeces
%
%> @retval A   Struct with Global matrices
%>
%==========================================================================

function [A] = LocalToGlobalMatrixElaElem(A, El, index)

A.M1_P_rho(index,index) = El.M1_P_rho_loc;
A.MPrjP_1(index,index)  = El.MPrjP_1_loc;

A.C1(index,index)       = El.C1_loc;
A.D1(index,index)       = El.D1_loc;

A.V1(index,index)       = El.V1_loc;
A.V2(index,index)       = El.V2_loc;
A.V3(index,index)       = El.V3_loc;
A.V4(index,index)       = El.V4_loc;
