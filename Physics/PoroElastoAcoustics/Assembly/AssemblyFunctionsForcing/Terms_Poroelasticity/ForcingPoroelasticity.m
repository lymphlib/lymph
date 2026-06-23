%> @file   ForcingPoroelasticity.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief  Global forcing structure for poroelasticity equation with SIP implementation.
%>
%==========================================================================
%> @section classForcingPoroelasticity Class description
%==========================================================================
%> @brief Global forcing structure for poroelasticity equation with SIP implementation.
%>
%> @param Forcing_loc   Forcing term struct from assembly function
%>
%> @retval Forcing      Final forcing term for assembly
%>                   
%==========================================================================

function [Forcing] = ForcingPoroelasticity(Forcing_loc)

    %% Symmetric Interior Penalty Method Forcing Term
    Forcing.f_p = [Forcing_loc.Volume.F1_loc; Forcing_loc.Volume.F2_loc];
    Forcing.g_p = [Forcing_loc.Volume.G1_loc; Forcing_loc.Volume.G2_loc];
    
    Forcing.j_p = [Forcing_loc.Volume.J1_loc; Forcing_loc.Volume.J2_loc];
    Forcing.h_p = [Forcing_loc.Volume.H1_loc; Forcing_loc.Volume.H2_loc];
    
    Forcing.f_p_diri = [Forcing_loc.Faces.F1_D_loc; Forcing_loc.Faces.F2_D_loc];
    Forcing.g_p_diri = [Forcing_loc.Faces.G1_D_loc; Forcing_loc.Faces.G2_D_loc];

end
