%> @file   VolumeForcingAssemblyPoroElastoAcousticsST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief Assembly of the forcing terms associated with volume integrals for
%> poroelastoacoustics problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyPoroElastoAcousticsST description
%==========================================================================
%> @brief            Assembly of the forcing terms associated with volume 
%> integrals for poroelastoacoustics problem (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Forcing   Forcing struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Forcing  Forcing struct
%>                   
%==========================================================================

function [Forcing] = VolumeForcingAssemblyPoroElastoAcousticsST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
      
  switch AssembInfo.label(ie)
      case 'P'
          %% Assembly of volume forcing terms for the poroelasticity
          [Forcing.Poro] = VolumeForcingAssemblyPoroelasticityST(Data, Forcing.Poro, elem, ie, id, nbases, AssembInfo);
      case 'E'
          %% Assembly of volume forcing terms for the elastodynamics
          [Forcing.Ela] = VolumeForcingAssemblyElastodynamicsST_PEA(Data, Forcing.Ela, elem, ie, id, nbases, AssembInfo);
      case 'A'
          %% Assembly of volume forcing terms for the acoustics
          [Forcing.Acu] = VolumeForcingAssemblyAcousticsST(Data, Forcing.Acu, elem, ie, id, nbases, AssembInfo);
  end

end
                    
