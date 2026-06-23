%> @file   ForcingLaplacian.m
%> @author Mattia Corti
%> @date   29 January 2026
%> @brief  Global forcing structure for Laplacian problem with SIP implementation.
%>
%==========================================================================
%> @section classForcingLaplacian Class description
%==========================================================================
%> @brief Global forcing structure for Laplacian problem with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing            Final forcing term for assembly
%>                   
%==========================================================================

function [Forcing] = ForcingLaplacian(Forcing_loc)

    %% Symmetric Interior Penalty Method Forcing Term
    Forcing.F = Forcing_loc.Volume.F_loc + Forcing_loc.Faces.F_D_loc;

    Forcing.Adaptivity = Forcing_loc;

end
