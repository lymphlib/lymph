%> @file   VolumeForcingAssemblyElastodynamicsST_PEA.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief Assembly of the vectors associated with volume integrals for
%> Elastodynamics equation (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyElastodynamicsST_PEA Class description
%==========================================================================
%> @brief            Assembly of the vectors associated with volume
%> integrals for Elastodynamics equation (subtriangulation implementation).
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

function [Forcing] = VolumeForcingAssemblyElastodynamicsST_PEA(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        fSource1 = Data.source_ue{1}(elem.xq,elem.yq);
        fSource2 = Data.source_ue{2}(elem.xq,elem.yq);
        gSource1 = Data.source_ued{1}(elem.xq,elem.yq);
        gSource2 = Data.source_ued{2}(elem.xq,elem.yq);
        
        MxxSource = Data.sourceMxx_el{1}(elem.xq,elem.yq);
        MxySource = Data.sourceMxy_el{1}(elem.xq,elem.yq);
        MyxSource = Data.sourceMyx_el{1}(elem.xq,elem.yq);
        MyySource = Data.sourceMyy_el{1}(elem.xq,elem.yq);

        % Forcing term assembly
        Forcing.F1_loc(1:nbases,1) = (elem.dx .* elem.phiq)'*fSource1 ...
                                       + (elem.dx .* elem.gradqx)'*MxxSource ...
                                       + 0.5 * (elem.dx .* elem.gradqy)'*MxySource;

        Forcing.F2_loc(1:nbases,1) = (elem.dx .* elem.phiq)'*fSource2 ...
                                       + 0.5 * (elem.dx .* elem.gradqx)'*MyxSource ...
                                       + (elem.dx .* elem.gradqy)'*MyySource;

        Forcing.G1_loc(1:nbases,1) = (elem.dx .* elem.phiq)'*gSource1;
        Forcing.G2_loc(1:nbases,1) = (elem.dx .* elem.phiq)'*gSource2;
        
end

