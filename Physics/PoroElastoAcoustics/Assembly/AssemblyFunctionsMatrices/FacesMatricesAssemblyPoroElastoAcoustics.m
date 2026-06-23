%> @file   FacesMatricesAssemblyPoroElastoAcoustics.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Assembly of the matrices associated with face's integrals for
%> poroelastoacoustics problem.
%>
%==========================================================================
%> @section classFacesMatricesAssemblyPoroElastoAcoustics Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for poroelastoacoustics problem.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Matrices struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = FacesMatricesAssemblyPoroElastoAcoustics(Data, femregion, Matrices, face, AssembInfo)
   
    switch AssembInfo.label(face.ie)
        case 'P'
            %% Assembly of face matrices for the poroelasticity
            [Matrices.Poro] = FacesMatricesAssemblyPoroelasticity(Data, femregion, Matrices.Poro, face, AssembInfo);
        
            if face.neigh_ie > 0
                %% Assembly of face matrices for the poroelasticity-elasticity coupling
                if AssembInfo.label(face.neigh_ie) == 'E'
                    [Matrices.PoroEla] = FacesMatricesAssemblyPoroElastic(Data, femregion, Matrices.PoroEla, face, AssembInfo);
                
                %% Assembly of face matrices for the poroelasticity-acoustic coupling
                elseif AssembInfo.label(face.neigh_ie) == 'A'
                    [Matrices.PoroAcu] = FacesMatricesAssemblyPoroAcoustic(Data, femregion, Matrices.PoroAcu, face, AssembInfo);
                end

            end
        
        case 'E' 
            %% Assembly of face matrices for the elastodynamics
            [Matrices.Ela] = FacesMatricesAssemblyElastodynamics_PEA(Data, femregion, Matrices.Ela, face, AssembInfo);
                            
            %% Assembly of face matrices for the poroelasticity-acoustic coupling
            if face.neigh_ie > 0
                if AssembInfo.label(face.neigh_ie) == 'A'
                    [Matrices.ElaAcu] = FacesMatricesAssemblyPoroAcoustic(Data, femregion, Matrices.ElaAcu, face, AssembInfo);
                end
            end
        case 'A'
            %% Assembly of face matrices for the acoustics
            [Matrices.Acu] = FacesMatricesAssemblyAcoustics(Data, femregion, Matrices.Acu, face, AssembInfo);
    end

end