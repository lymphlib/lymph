%> @file   VolumeMatricesAssemblyPoroElastoAcousticsQF.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Assembly of the matrices associated with volume integrals for
%> poroelastoacoustics problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeMatricesAssemblyPoroElastoAcousticsQF description
%==========================================================================
%> @brief            Assembly of the matrices associated with volume 
%> integrals for poroelastoacoustics problem (QF implementation).
%>
%> @param Data       Struct with problem's data
%> @param Integral   Struct containing the integrals computed with the QF implementation
%> @param Matrices   Matrices struct (to be modified)
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = VolumeMatricesAssemblyPoroElastoAcousticsQF(Data, Integral, Matrices, ie, id, nbases, AssembInfo)
      
  switch AssembInfo.label(ie)
      case 'P'
          %% Assembly of volume matrices for the poroelasticity
          [Matrices.Poro] = VolumeMatricesAssemblyPoroelasticityQF(Data, Integral, Matrices.Poro, ie, id, nbases, AssembInfo);
      case 'E'
          %% Assembly of volume matrices for the elastodynamics (from Physics/Elastodynamics)
          [Matrices.Ela] = VolumeMatricesAssemblyElastodynamicsQF(Data, Integral, Matrices.Ela, ie, id, nbases, AssembInfo);
    case 'A'
          %% Assembly of volume matrices for the acoustics
          [Matrices.Acu] = VolumeMatricesAssemblyAcousticsQF(Data, Integral, Matrices.Acu, ie, id, nbases, AssembInfo);
  end

end
                    
