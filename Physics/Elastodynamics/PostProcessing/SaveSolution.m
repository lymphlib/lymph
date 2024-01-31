%> @file  SaveSolution.m
%> @author Ilario Mazzieri
%> @date 11 May 2023
%> @brief  Save solution vector in a struct
%>
%==========================================================================
%> @section classSaveSolution Class description
%> @brief  Save solution vector in a struct
%
%> @param Uh     Solution vector
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Solutions  Struct with solution vector
%>
%==========================================================================

function  [Solutions]  = SaveSolution(Uh,femregion)

ne = femregion.ndof_e;


% elastic
ue_h   = Uh(1 : 2*ne);

ndof_vel = 2*ne;
dot_ue_h   = Uh(ndof_vel+1  : ndof_vel+2*ne);

Solutions.ue_h      = ue_h;
Solutions.dot_ue_h  = dot_ue_h;

