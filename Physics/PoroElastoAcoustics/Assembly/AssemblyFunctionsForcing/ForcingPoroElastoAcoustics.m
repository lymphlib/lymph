%> @file   ForcingPoroElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief  Global forcing terms structure for poroelastoacoustics problem with SIP implementation.
%>
%==========================================================================
%> @section classForcingPoroElastoAcoustics Class description
%==========================================================================
%> @brief           Global forcing terms structure for poroelastoacoustics problem with SIP implementation.
%>
%> @param Forcing_loc  Forcing struct from assembly function
%>
%> @retval Forcing     Forcing struct containing global forcing terms for solver
%>                   
%==========================================================================

function [Forcing] = ForcingPoroElastoAcoustics(Forcing_loc)

    %% SIP forcing terms for the elastodynamics
    Forcing_Ela_loc.Volume = Forcing_loc.Volume.Ela;
    Forcing_Ela_loc.Faces  = Forcing_loc.Faces.Ela;
    [Forcing.Ela] = ForcingElastodynamics_PEA(Forcing_Ela_loc);
    
    %% SIP forcing terms for the poroelasticity
    Forcing_Poro_loc.Volume = Forcing_loc.Volume.Poro;
    Forcing_Poro_loc.Faces  = Forcing_loc.Faces.Poro;
    [Forcing.Poro] = ForcingPoroelasticity(Forcing_Poro_loc);

    %% SIP forcing terms for the acoustics
    Forcing_Acu_loc.Volume = Forcing_loc.Volume.Acu;
    Forcing_Acu_loc.Faces  = Forcing_loc.Faces.Acu;
    [Forcing.Acu] = ForcingAcoustics(Forcing_Acu_loc);

end
