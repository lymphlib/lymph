%> @file   SIPMatricesFK.m
%> @author Mattia Corti
%> @date   16 September 2025
%> @brief  Global matrices structure for FK problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesFK Class description
%==========================================================================
%> @brief           Global matrices structure for FK problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices  Matrices struct containing global sparse matrices
%>                   
%==========================================================================

function [Matrices] = SIPMatricesFK(Matrices_loc)

    %% Symmetric Interior Penalty Method Matrices

    % Projection matrix
    Matrices.M_prj = Matrices_loc.Volume.M_prj_loc;
    
    % Mass matrix associated with the reaction term
    Matrices.M    = Matrices_loc.Volume.M_loc;

    % IPDG stiffness matrix
    Matrices.dGA  = Matrices_loc.Volume.A_loc + Matrices_loc.Faces.SA_loc;
    Matrices.A    = Matrices.dGA - Matrices_loc.Faces.IA_loc - transpose(Matrices_loc.Faces.IA_loc);   % DG Stiffness Matrix

    % Matrix associated with the trilinear term
    Matrices.NonLinear.M_NL = Matrices_loc.Volume3L.M_NL_loc;

    %% Save matrices for adaptivity
    Matrices.Adaptivity = Matrices_loc;
end
