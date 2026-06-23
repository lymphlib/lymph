%> @file   SIPMatricesElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   28 May 2026
%> @brief  Global matrices structure for elastodynamics-acoustic coupling
%> terms.
%>
%==========================================================================
%> @section classSIPMatricesElastoAcoustics Class description
%==========================================================================
%> @brief           Global matrices structure for elastodynamics-acoustic coupling terms.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesElastoAcoustics(Matrices_loc)

    %% Matrices
    Matrices.C1_E = [Matrices_loc.Faces.C1_loc ; Matrices_loc.Faces.C2_loc];
    Matrices.C2_E = [Matrices_loc.Faces.C1_loc ; Matrices_loc.Faces.C2_loc];

    Matrices.C1_A = - Matrices.C1_E';
    Matrices.C2_A = - Matrices.C2_E';

end