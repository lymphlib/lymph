%> @file  LocalToGlobalMatrixElaEdge.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief Assembly of the edge matrices for the elastic domain
%>
%==========================================================================
%> @section classLocalToGlobalMatrixElaEdge Class description
%==========================================================================
%> @brief Assembly of the edge matrices for the elastic domain
%
%> @param A      Struct with global matrices
%> @param Edge   Struct with local edge matrices
%> @param index  Vector wuth local to global mapping
%
%> @retval A   Struct with global matrices
%>
%==========================================================================

function  [A] = LocalToGlobalMatrixElaEdge(A, Edge, index)

A.S1_P(index,index)   = Edge.S1_P_loc;
A.S4_P(index,index)   = Edge.S4_P_loc;

A.IT1_P(index,index)  = Edge.IT1_P_loc;
A.IT2_P(index,index)  = Edge.IT2_P_loc;
A.IT3_P(index,index)  = Edge.IT3_P_loc;
A.IT4_P(index,index)  = Edge.IT4_P_loc;

A.ABC_R1(index,index) = Edge.ABC_R1_loc;
A.ABC_R2(index,index) = Edge.ABC_R2_loc;
A.ABC_R3(index,index) = Edge.ABC_R3_loc;
A.ABC_R4(index,index) = Edge.ABC_R4_loc;

A.ABC_S1(index,index) = Edge.ABC_S1_loc;
A.ABC_S2(index,index) = Edge.ABC_S2_loc;
A.ABC_S3(index,index) = Edge.ABC_S3_loc;
A.ABC_S4(index,index) = Edge.ABC_S4_loc;
