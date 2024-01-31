%> @file  AssembleODEMatrices.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Assembly of the matrices A,B and C for the problem Ax'' + Bx' + Cx = F
%>
%==========================================================================
%> @section classAssembleODEMatrices Class description
%> @brief  Assembly of the matrices A,B and C for the problem Ax'' + Bx' + Cx = F
%
%> @param Matrices   Matrices.Ela     = Matrices for elastic domain
%
%> @retval A,B,C     Matrices for the problem Ax'' + Bx' + Cx = F
%>
%==========================================================================

function [A, B, C] = AssembleODEMatrices(Matrices)  
                            

A = Matrices.Ela.M_P_rho;
 
B = Matrices.Ela.Dvel + Matrices.Ela.Svel;

C = Matrices.Ela.A_E + Matrices.Ela.Ddis + Matrices.Ela.Rdis;
        
 end

