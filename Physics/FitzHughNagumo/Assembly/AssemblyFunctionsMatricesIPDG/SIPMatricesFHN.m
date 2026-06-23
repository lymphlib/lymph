%> @file   SIPMatricesFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief  Global matrices structure for FHN problem with SIP implementation.
%>
%==========================================================================
%> @section classSIPMatricesFHN Class description
%==========================================================================
%> @brief           Global matrices structure for FHN problem with SIP implementation.
%>
%> @param Matrices_loc  Matrices struct from assembly function
%>
%> @retval Matrices  Matrices struct containing global sparse matrices
%>                   
%==========================================================================

function [Matrices] = SIPMatricesFHN(Matrices_loc)

    %% Symmetric Interior Penalty Method Matrices

    % Projection matrix
    Matrices.M_prj = Matrices_loc.Volume.M_prj_loc;
    
    % Mass matrix associated with the reaction term
    Matrices.M_u    = Matrices_loc.Volume.M_u_loc;
    Matrices.M_w    = Matrices_loc.Volume.M_w_loc;

    % IPDG stiffness matrix
    Matrices.dGA  = Matrices_loc.Volume.A_loc + Matrices_loc.Faces.SA_loc;
    Matrices.A    = Matrices.dGA - Matrices_loc.Faces.IA_loc - transpose(Matrices_loc.Faces.IA_loc);   % DG Stiffness Matrix

    %% Save matrices for adaptivity
    Matrices.Adaptivity = Matrices_loc;
end
