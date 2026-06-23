%> @file   VolumeForcingAssemblyPoroelasticityST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief Assembly of the vectors associated with volume integrals for
%> poroelasticity equation (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyPoroelasticityST Class description
%==========================================================================
%> @brief            Assembly of the vectors associated with volume
%> integrals for poroelasticity equation (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Forcing    Forcing struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = VolumeForcingAssemblyPoroelasticityST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters            
        fSource1 = Data.source_up{1}(elem.xq,elem.yq);
        fSource2 = Data.source_up{2}(elem.xq,elem.yq);
        jSource1 = Data.source_upd{1}(elem.xq,elem.yq);
        jSource2 = Data.source_upd{2}(elem.xq,elem.yq);

        gSource1 = Data.source_wp{1}(elem.xq,elem.yq);
        gSource2 = Data.source_wp{2}(elem.xq,elem.yq);
        hSource1 = Data.source_wpd{1}(elem.xq,elem.yq);
        hSource2 = Data.source_wpd{2}(elem.xq,elem.yq);

        MxxSource = Data.sourceMxx_poro{1}(elem.xq,elem.yq);
        MxySource = Data.sourceMxy_poro{1}(elem.xq,elem.yq);
        MyxSource = Data.sourceMyx_poro{1}(elem.xq,elem.yq);
        MyySource = Data.sourceMyy_poro{1}(elem.xq,elem.yq);

        % Vector assembling
        Forcing.F1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * fSource1 ...
                                   + (elem.dx .* elem.gradqx)' * MxxSource ...
                                   + 0.5 * (elem.dx .* elem.gradqy)' * MxySource;
        Forcing.F2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * fSource2 ...
                                   + 0.5 * (elem.dx .* elem.gradqx)' * MyxSource ...
                                   + (elem.dx .* elem.gradqy)' * MyySource;

        Forcing.J1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * jSource1;
        Forcing.J2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * jSource2;

        Forcing.G1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * gSource1 ...
                                   + (elem.dx .* elem.gradqx)'*MxxSource ...
                                   + 0.5 * (elem.dx .* elem.gradqy)'*MxySource;
        Forcing.G2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * gSource2 ...
                                   + 0.5 * (elem.dx .* elem.gradqx)' * MyxSource ...
                                   + (elem.dx .* elem.gradqy)' * MyySource;
        
        Forcing.H1_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * hSource1;
        Forcing.H2_loc(1:nbases,1) = (elem.dx .* elem.phiq)' * hSource2;

end

