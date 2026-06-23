%> @file   SIPMatricesHeat.m
%> @author Mattia Corti
%> @date   5 February 2026
%> @brief  Global matrices structure for heat equation with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesHeat Class description
%==========================================================================
%> @brief           Global matrices structure for heat equation with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesHeat(Matrices_loc)

    %% Symmetric Interior Penalty Method Matrices
    Matrices.Mprj = Matrices_loc.Volume.Mprj_loc;
    Matrices.M    = Matrices_loc.Volume.M_loc;
    
    Matrices.dGA  = Matrices_loc.Volume.A_loc + Matrices_loc.Faces.SA_loc;
    Matrices.A    = Matrices.dGA - Matrices_loc.Faces.IA_loc - transpose(Matrices_loc.Faces.IA_loc);   % DG Stiffness Matrix

    %% Save matrices for adaptivity
    Matrices.Adaptivity = Matrices_loc;
end