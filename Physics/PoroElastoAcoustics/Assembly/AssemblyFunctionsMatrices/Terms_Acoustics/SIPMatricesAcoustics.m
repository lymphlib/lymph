%> @file   SIPMatricesAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   12 May 2026
%> @brief  Global matrices structure for acoustics problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesAcoustics Class description
%==========================================================================
%> @brief           Global matrices structure for acoustics problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesAcoustics(Matrices_loc)

    %% Mass matrices
    Matrices.MPrjA = Matrices_loc.Volume.MPrjA_loc;
    Matrices.M_A   = Matrices_loc.Volume.M_P_A_loc;

    %% Stiffness matrices
    Matrices.A_A = Matrices_loc.Volume.W_loc + Matrices_loc.Faces.S_A_loc - Matrices_loc.Faces.IT_A_loc - Matrices_loc.Faces.IT_A_loc';
    Matrices.dGA = Matrices_loc.Volume.W_loc + Matrices_loc.Faces.S_A_loc;

end