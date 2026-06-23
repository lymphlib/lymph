%> @file   ForcingElastodynamics_PEA.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief  Global forcing structure for elastodynamics equation with SIP implementation.
%>
%==========================================================================
%> @section classForcingElastodynamics_PEA Class description
%==========================================================================
%> @brief Global forcing structure for elastodynamics equation with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing      Final forcing term for assembly
%>                   
%==========================================================================

function [Forcing] = ForcingElastodynamics_PEA(Forcing_loc)

    %% Symmetric Interior Penalty Method Forcing Term
    Forcing.f_e = [Forcing_loc.Volume.F1_loc; Forcing_loc.Volume.F2_loc];
    Forcing.f_e_diri = [Forcing_loc.Faces.F1_D_loc; Forcing_loc.Faces.F2_D_loc];

    Forcing.g_e = [Forcing_loc.Volume.G1_loc; Forcing_loc.Volume.G2_loc];

end
