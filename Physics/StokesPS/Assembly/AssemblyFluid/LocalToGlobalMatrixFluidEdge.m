%> @file  LocalToGlobalMatrixFluidEdge.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Assembly of the edge matrices for the fluid domain
%>
%==========================================================================
%> @section classLocalToGlobalMatrixFluidEdge Class description
%==========================================================================
%> @brief Assembly of the edge matrices for the fluid domain
%>
%> @param A      Struct with global matrices
%> @param Edge   Struct with local edge matrices
%> @param index  Vector wuth local to global mapping
%> 
%> @retval A   Struct with global matrices
%>
%==========================================================================

function  [A] = LocalToGlobalMatrixFluidEdge(A, Edge, index)

A.S1_B(index,index)          = Edge.S1_B_loc;
A.S2_B(index,index)          = Edge.S2_B_loc;
A.S3_B(index,index)          = Edge.S3_B_loc;
A.S4_B(index,index)          = Edge.S4_B_loc;
A.BT1(index,index)           = Edge.BT1_loc;
A.BT2(index,index)           = Edge.BT2_loc;
A.BT3(index,index)           = Edge.BT3_loc;
A.BT4(index,index)           = Edge.BT4_loc;
A.C1(index,index)            = Edge.C1_loc;
A.C2(index,index)            = Edge.C2_loc;
A.C3(index,index)            = Edge.C3_loc;
A.C4(index,index)            = Edge.C4_loc;
A.C5(index,index)            = Edge.C5_loc;
A.C6(index,index)            = Edge.C6_loc;
A.C7(index,index)            = Edge.C7_loc;
A.C8(index,index)            = Edge.C8_loc;
A.C9(index,index)            = Edge.C9_loc;
