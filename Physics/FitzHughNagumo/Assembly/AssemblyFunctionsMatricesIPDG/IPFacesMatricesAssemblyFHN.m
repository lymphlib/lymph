%> @file  IPFacesMatricesAssemblyFHN.m
%> @author Mattia Corti, Cateerina Leimer Saglio
%> @date   11 May 2026
%> @brief Assembly of the matrices associated with face integrals for FHN
%> problem with IPDG dicretization.
%>
%==========================================================================
%> @section classIPFacesMatricesAssemblyFHN Class description
%==========================================================================
%> @brief            Assembly of the matrices associated with face 
%> integrals for FHN problem with IPDG dicretization.
%>
%> @param Data       Struct with problem's data
%> @param femregion  Finite Element struct (see CreateDOF.m)
%> @param Matrices   Matrices struct (to be modified)
%> @param face       Face struct containing the quadrature nodes, penalty and bases
%> @param AssembInfo Struct containing specific information for the assembly procedure
%>
%> @retval Matrices  Matrices struct
%>                   
%==========================================================================

function [Matrices] = IPFacesMatricesAssemblyFHN(Data, femregion, Matrices, face, AssembInfo)
   
        ie = face.ie;
        iedg = face.iedg;
        bas_ie = 1:femregion.nbases(ie);

        % Evaluation of physical parameters
        D_ext  = Data.D_ext{femregion.id(ie)}(face.xq,face.yq,0);
        
        if ~Data.isotropy
            D_axn  = Data.D_axn{femregion.id(ie)}(face.xq,face.yq,0);

            % Axonal directions
            D_axn_xx = AssembInfo.Tens.D_xx(ie);
            D_axn_xy = AssembInfo.Tens.D_xy(ie);
            D_axn_yx = AssembInfo.Tens.D_yx(ie);
            D_axn_yy = AssembInfo.Tens.D_yy(ie);
        end
    
    %% Dirichlet boundary faces
    if face.neigh_ie == -1

        if ~Data.isotropy
            penalty = (D_ext+D_axn) .* face.penalty_geom.harm(iedg);
        else
            penalty = D_ext .* face.penalty_geom.harm(iedg);
        end

        % Extracellular component
        Matrices.IA_loc(bas_ie,bas_ie) = Matrices.IA_loc(bas_ie,bas_ie) ...
                                            + (face.ds .* D_ext .* (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy))' * face.phiedgeq;
        Matrices.SA_loc(bas_ie,bas_ie) = Matrices.SA_loc(bas_ie,bas_ie) ...
                                            + (face.ds .* penalty .* face.phiedgeq)' * face.phiedgeq;

        % Axonal Component
        if ~Data.isotropy
            Matrices.IA_loc(bas_ie,bas_ie) = Matrices.IA_loc(bas_ie,bas_ie) ...
                                                + (face.ds .* ((D_axn .* (D_axn_xx * face.nx + D_axn_yx * face.ny)) .* face.gradedgeqx))' * face.phiedgeq ...
                                                + (face.ds .* ((D_axn .* (D_axn_xy * face.nx + D_axn_yy * face.ny)) .* face.gradedgeqy))' * face.phiedgeq;
        end
    
        %% Internal faces
    elseif face.neigh_ie > 0

        % Evaluation of physical parameters
        D_ext_n  = Data.D_ext{femregion.id(face.neigh_ie)}(face.xq,face.yq,0);
        
        if ~Data.isotropy
            D_axn_n = Data.D_axn{femregion.id(face.neigh_ie)}(face.xq,face.yq,0);
            penalty = 2./(1./(D_ext+D_axn)+1./(D_ext_n+D_axn_n)) .* face.penalty_geom.harm(iedg);
        else
            penalty = 2./(1./(D_ext)+1./(D_ext_n)) .* face.penalty_geom.harm(iedg);
        end

        %% Element itself

        % Extracellular component
        Matrices.IA_loc(bas_ie,bas_ie) = Matrices.IA_loc(bas_ie,bas_ie) ...
                                            + 0.5 * (face.ds .* D_ext .* (face.nx * face.gradedgeqx + face.ny * face.gradedgeqy))' * face.phiedgeq;
        Matrices.SA_loc(bas_ie,bas_ie) = Matrices.SA_loc(bas_ie,bas_ie) ...
                                            + (face.ds .* penalty .* face.phiedgeq)' * face.phiedgeq;

        % Axonal Component
        if ~Data.isotropy
            Matrices.IA_loc(bas_ie,bas_ie) = Matrices.IA_loc(bas_ie,bas_ie) ...
                                                + 0.5 * (face.ds .* ((D_axn .* (D_axn_xx * face.nx + D_axn_yx * face.ny)) .* face.gradedgeqx))' * face.phiedgeq ...
                                                + 0.5 * (face.ds .* ((D_axn .* (D_axn_xy * face.nx + D_axn_yy * face.ny)) .* face.gradedgeqy))' * face.phiedgeq;
        end


        %% Neighboring element
        bas_neigh_ie = 1:femregion.nbases(face.neigh_ie);

        % Extracellular component
        Matrices.IA_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IA_loc(face.neigh_idx,bas_neigh_ie) ...
                                                          - 0.5 * (face.ds .* D_ext .* ( face.nx * face.gradedgeqx +  face.ny * face.gradedgeqy))' * face.phiedgeqneigh;
        
        Matrices.SA_loc(face.neigh_idx,bas_neigh_ie) = Matrices.SA_loc(face.neigh_idx,bas_neigh_ie) ...
                                                          - (face.ds .* penalty .* face.phiedgeq)' * face.phiedgeqneigh;

        % Axonal Component
        if ~Data.isotropy
            Matrices.IA_loc(face.neigh_idx,bas_neigh_ie) = Matrices.IA_loc(face.neigh_idx,bas_neigh_ie) ...
                                                                - 0.5 * (face.ds .* ((D_axn .* (D_axn_xx * face.nx + D_axn_yx * face.ny)) .* face.gradedgeqx))' * face.phiedgeqneigh ...
                                                                - 0.5 * (face.ds .* ((D_axn .* (D_axn_xy * face.nx + D_axn_yy * face.ny)) .* face.gradedgeqy))' * face.phiedgeqneigh;
        end

    end

end
