%> @file   ForcingPreallocationPoroElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Preallocation of the forcing terms structure for the poroelastoacoustics problem.
%>
%==========================================================================
%> @section classForcingPreallocationPoroElastoAcoustics Class description
%==========================================================================
%> @brief           Preallocation of the forcing terms structure for poroelastoacoustics 
%> problem.
%>
%> @param GenForcing Struct containing the possible forcing terms to assemble:
%>                     - Vector: Preallocated vector
%>
%> @retval Forcing  Forcing struct containing local forcing terms stored using
%> cells associated with the mesh elements
%>                   
%==========================================================================

function [Forcing] = ForcingPreallocationPoroElastoAcoustics(GenForcing)

    %% Preallocation of the forcing terms for the elastodynamics
    [ForcingEla] = ForcingPreallocationElastodynamics(GenForcing);
    
    %% Preallocation of the forcing terms for the poroelasticity
    [ForcingPoro] = ForcingPreallocationPoroelasticity(GenForcing);

    %% Preallocation of the forcing terms for the acoustics
    [ForcingAcu] = ForcingPreallocationAcoustics(GenForcing);

    %% Creation of final forcing terms structure
    Forcing.Volume.Ela  = ForcingEla.Volume;
    Forcing.Volume.Poro = ForcingPoro.Volume;
    Forcing.Volume.Acu  = ForcingAcu.Volume;

    Forcing.Faces.Ela  = ForcingEla.Faces;
    Forcing.Faces.Poro = ForcingPoro.Faces;
    Forcing.Faces.Acu  = ForcingAcu.Faces;

end
    