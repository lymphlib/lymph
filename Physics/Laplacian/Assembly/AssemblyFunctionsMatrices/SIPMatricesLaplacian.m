%> @file   SIPMatricesLaplacian.m
%> @author Mattia Corti
%> @date   29 January 2026
%> @brief  Global matrices structure for laplacian problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesLaplacian Class description
%==========================================================================
%> @brief           Global matrices structure for laplacian problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices     Matrices struct containing global matrices for solver
%>                   
%==========================================================================

function [Matrices] = SIPMatricesLaplacian(Matrices_loc)

    %% Symmetric Interior Penalty Method Matrices
    Matrices.Mprj = Matrices_loc.Volume.Mprj_loc;
    
    Matrices.dGA  = Matrices_loc.Volume.A_loc + Matrices_loc.Faces.SA_loc;
    Matrices.A    = Matrices.dGA - Matrices_loc.Faces.IA_loc - transpose(Matrices_loc.Faces.IA_loc);   % DG Stiffness Matrix

    %% Save matrices for adaptivity
    Matrices.Adaptivity = Matrices_loc;
end