%> @file   FacesForcingAssemblyPoroElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief Assembly of the forcing terms associated with face's integrals for
%> poroelastoacoustics problem.
%>
%==========================================================================
%> @section classFacesForcingAssemblyPoroElastoAcoustics Class description
%==========================================================================
%> @brief            Assembly of the forcing terms associated with face's 
%> integrals for poroelastoacoustics problem.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Forcing   Forcing struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Forcing  Forcing struct
%>                   
%==========================================================================

function [Forcing] = FacesForcingAssemblyPoroElastoAcoustics(Data, femregion, Forcing, face, AssembInfo)
   
    switch AssembInfo.label(face.ie)
        case 'P'
            %% Assembly of face forcing terms for the poroelasticity
            [Forcing.Poro] = FacesForcingAssemblyPoroelasticity(Data, femregion, Forcing.Poro, face, AssembInfo);
        
        case 'E' 
            %% Assembly of face forcing terms for the elastodynamics
            [Forcing.Ela] = FacesForcingAssemblyElastodynamics_PEA(Data, femregion, Forcing.Ela, face, AssembInfo);
                            
        case 'A'
            %% Assembly of face forcing terms for the acoustics
            [Forcing.Acu] = FacesForcingAssemblyAcoustics(Data, femregion, Forcing.Acu, face, AssembInfo);
    end

end