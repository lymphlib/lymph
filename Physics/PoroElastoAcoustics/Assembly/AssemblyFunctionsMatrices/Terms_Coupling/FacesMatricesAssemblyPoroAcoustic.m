%> @file   FacesMatricesAssemblyPoroAcoustic.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   29 May 2026
%> @brief Assembly of the matrices associated with face's integrals for
%> poroelastic-acoustic coupling.
%>
%==========================================================================
%> @section classFacesMatricesAssemblyPoroAcoustic Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face's 
%> integrals for poroelastic-acoustic coupling.
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

function [Matrices] = FacesMatricesAssemblyPoroAcoustic(Data, femregion, Matrices, face, AssembInfo)

    % Evaluation of physical parameters
    rho_a = Data.rho_a{femregion.id_phys(face.neigh_ie)}(face.xq,face.yq);
    
    bas_neigh_ie = 1:femregion.nbases(face.neigh_ie);
        
    Matrices.C1_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C1_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    + (face.ds .* rho_a .* face.nx .* face.phiedgeq)' * face.phiedgeqneigh;

    Matrices.C2_loc(face.neigh_idx,bas_neigh_ie) = Matrices.C2_loc(face.neigh_idx,bas_neigh_ie) ...
                                                    + (face.ds .* rho_a .* face.ny .* face.phiedgeq)' * face.phiedgeqneigh;
                                                 
end