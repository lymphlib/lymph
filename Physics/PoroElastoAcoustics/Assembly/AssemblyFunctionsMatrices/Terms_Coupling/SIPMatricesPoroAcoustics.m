%> @file   SIPMatricesPoroAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief  Global matrices structure for poroelasticity-acoustic coupling.
%>
%==========================================================================
%> @section classSIPMatricesPoroAcoustics Class description
%==========================================================================
%> @brief           Global matrices structure for poroelasticity-acoustic coupling.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesPoroAcoustics(Matrices_loc)

    %% Matrices
    Matrices.C1_P = [Matrices_loc.Faces.C1_loc ; Matrices_loc.Faces.C2_loc];
    Matrices.C2_P = [Matrices_loc.Faces.C1_loc ; Matrices_loc.Faces.C2_loc];

    Matrices.C1_A = - Matrices.C1_P';
    Matrices.C2_A = - Matrices.C2_P';

end