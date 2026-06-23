%> @file   VolumeMatricesAssemblyPoroElastoAcousticsST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> poroelastoacoustics problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyPoroElastoAcousticsST description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for poroelastoacoustics problem (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Matrices   Matrices struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = VolumeMatricesAssemblyPoroElastoAcousticsST(Data, Matrices, elem, ie, id, nbases, AssembInfo)
      
  switch AssembInfo.label(ie)
      case 'P'
          %% Assembly of volume matrices for the poroelasticity
          [Matrices.Poro] = VolumeMatricesAssemblyPoroelasticityST(Data, Matrices.Poro, elem, ie, id, nbases, AssembInfo);
      case 'E'
          %% Assembly of volume matrices for the elastodynamics (from Physics/Elastodynamics)
          [Matrices.Ela] = VolumeMatricesAssemblyElastodynamicsST(Data, Matrices.Ela, elem, ie, id, nbases, AssembInfo);
      case 'A'
          %% Assembly of volume matrices for the acoustics
          [Matrices.Acu] = VolumeMatricesAssemblyAcousticsST(Data, Matrices.Acu, elem, ie, id, nbases, AssembInfo);
  end

end
                    
