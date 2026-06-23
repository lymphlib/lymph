%> @file   FacesForcingAssemblyElastodynamics.m
%> @author Mattia Corti
%> @date   12 May 2026
%> @brief Assembly of the forcing vectors associated with face's integrals for
%> elastodynamics equation.
%>
%==========================================================================
%> @section classFacesForcingAssemblyElastodynamics Class description
%==========================================================================
%> @brief            Assembly of the forcing vectors associated with face's
%> integrals for elastodynamics equation.
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

function [Forcing] = FacesForcingAssemblyElastodynamics(Data, femregion, Forcing, face, AssembInfo)
   
    ie   = face.ie;
    iedg = face.iedg;

    % Evaluation of physical parameters
    mu  = Data.mu_el{femregion.id(ie)}(face.xq,face.yq);
    lam = Data.lam_el{femregion.id(ie)}(face.xq,face.yq);

    % Evaluation of boundary functions
    gD1 = Data.DirBCEla{1}(face.xq,face.yq,AssembInfo.t);
    gD2 = Data.DirBCEla{2}(face.xq,face.yq,AssembInfo.t);

    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        Forcing.F1_D_loc(1:femregion.nbases(ie),1) = Forcing.F1_D_loc(1:femregion.nbases(ie),1) ...
                                                       - (face.ds.* ( ((lam+2*mu) * face.nx) .* face.gradedgeqx + (mu * face.ny) .* face.gradedgeqy))' * gD1 ...
                                                       - (face.ds.* ((lam * face.ny) .* face.gradedgeqx + (mu * face.nx) .* face.gradedgeqy))' * gD2 ...
                                                       + face.penalty_geom.harm(iedg) * (face.ds .* (lam + 2*mu) .* face.phiedgeq)' * gD1;

        Forcing.F2_D_loc(1:femregion.nbases(ie),1) = Forcing.F2_D_loc(1:femregion.nbases(ie),1) ...
                                                       - (face.ds.* ((mu * face.ny) .* face.gradedgeqx + (lam * face.nx) .* face.gradedgeqy))' * gD1 ...
                                                       - (face.ds.* ((mu * face.nx) .* face.gradedgeqx + ((lam+2*mu) * face.ny) .* face.gradedgeqy))' * gD2 ...
                                                       + face.penalty_geom.harm(iedg) * (face.ds .* (lam + 2*mu) .* face.phiedgeq)' * gD2;        
    end

end
