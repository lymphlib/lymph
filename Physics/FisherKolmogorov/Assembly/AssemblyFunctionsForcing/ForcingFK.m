%> @file   ForcingFK.m
%> @author Mattia Corti
%> @date   14 April 2026
%> @brief  Global forcing structure for FK problem with SIP implementation.
%>
%==========================================================================
%> @section classForcingFK Class description
%==========================================================================
%> @brief           Global forcing structure for FK problem with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing      Final forcing term for assembly
%>                   
%==========================================================================

function [Forcing] = ForcingFK(Forcing_loc)

    % Symmetric Interior Penalty Method Forcing Term
    Forcing.F = Forcing_loc.Volume.F_loc + Forcing_loc.Faces.F_D_loc;
    
    Forcing.Adaptivity = Forcing_loc;

end
