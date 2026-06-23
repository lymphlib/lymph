%> @file   ForcingElastodynamics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   8 May 2026
%> @brief  Global forcing structure for elastodynamics equation with SIP implementation.
%>
%==========================================================================
%> @section classForcingElastodynamics Class description
%==========================================================================
%> @brief Global forcing structure for elastodynamics equation with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing      Final forcing term for assembly
%>                   
%==========================================================================

function [Forcing] = ForcingElastodynamics(Forcing_loc)

    %% Symmetric Interior Penalty Method Forcing Term
    Forcing.F = [Forcing_loc.Volume.F1_loc + Forcing_loc.Faces.F1_D_loc;
                 Forcing_loc.Volume.F2_loc + Forcing_loc.Faces.F2_D_loc];

    Forcing.G = [Forcing_loc.Volume.G1_loc; Forcing_loc.Volume.G2_loc];

end
