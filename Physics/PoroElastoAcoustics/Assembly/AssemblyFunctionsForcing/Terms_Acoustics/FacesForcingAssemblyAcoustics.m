%> @file   FacesForcingAssemblyAcoustics.m
%> @author Mattia Corti
%> @date   30 May 2026
%> @brief Assembly of the forcing vectors associated with face's integrals for
%> acoustics equation.
%>
%==========================================================================
%> @section classFacesForcingAssemblyAcoustics Class description
%==========================================================================
%> @brief            Assembly of the forcing vectors associated with face's
%> integrals for acoustics equation.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Forcing    Forcing struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the
%> assembly procedure
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = FacesForcingAssemblyAcoustics(Data, femregion, Forcing, face, AssembInfo)
   

    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        ie   = face.ie;
        iedg = face.iedg;

        % Evaluation of physical parameters
        rho_a = Data.rho_a{femregion.id_phys(ie)}(face.xq,face.yq);
        gD1 = Data.DirBCAcu{1}(face.xq,face.yq);

        % Vector assembling
        Forcing.F1_D_loc(1:femregion.nbases(ie),1) = Forcing.F1_D_loc(1:femregion.nbases(ie),1) ...
                                                   - ( face.ds.* rho_a .* ( face.nx .* face.gradedgeqx + face.ny .* face.gradedgeqy))' * gD1 ...
                                                   + face.penalty_geom.harm(iedg) * (face.ds .* rho_a .* face.phiedgeq)' * gD1;

    end

end
