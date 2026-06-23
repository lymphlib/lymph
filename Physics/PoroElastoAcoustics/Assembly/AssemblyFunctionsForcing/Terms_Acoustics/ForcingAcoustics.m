%> @file   ForcingAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief  Global forcing structure for acoustics equation with SIP implementation.
%>
%==========================================================================
%> @section classForcingAcoustics Class description
%==========================================================================
%> @brief Global forcing structure for acoustics equation with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing      Final forcing term for assembly
%>                   
%==========================================================================

function [Forcing] = ForcingAcoustics(Forcing_loc)
    
    %% Symmetric Interior Penalty Method Forcing Term    
    Forcing.f_a = Forcing_loc.Volume.F1_loc;
    Forcing.f_a_diri = Forcing_loc.Faces.F1_D_loc;

end
