%> @file   ForcingFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date   11 May 2026
%> @brief  Global forcing structure for FHN problem with SIP implementation.
%>
%==========================================================================
%> @section classForcingFHN Class description
%==========================================================================
%> @brief           Global forcing structure for FHN problem with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing      Final forcing term for assembly
%>
%==========================================================================

function [Forcing] = ForcingFHN(Forcing_loc)

% Symmetric Interior Penalty Method Forcing Term
if length(Forcing_loc.Faces.F_D_loc) == length(Forcing_loc.Volume.F_loc)
    Forcing.Iext = Forcing_loc.Volume.Iext_loc + Forcing_loc.Faces.F_D_loc;
else
    Forcing.Iext = Forcing_loc.Volume.Iext_loc;
end

Forcing.F = Forcing_loc.Volume.F_loc;
Forcing.G = Forcing_loc.Volume.G_loc;


Forcing.Adaptivity = Forcing_loc;

end
