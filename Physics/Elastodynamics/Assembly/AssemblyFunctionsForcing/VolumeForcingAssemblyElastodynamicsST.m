%> @file   VolumeForcingAssemblyElastodynamicsST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   8 May 2026
%> @brief Assembly of the vectors associated with volume integrals for
%> Elastodynamics equation (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyElastodynamicsST Class description
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

function [Forcing] = VolumeForcingAssemblyElastodynamicsST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters            
        fSource1 = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t);
        fSource2 = Data.source_f{2}(elem.xq,elem.yq,AssembInfo.t);
        gSource1 = Data.source_g{1}(elem.xq,elem.yq,AssembInfo.t);
        gSource2 = Data.source_g{2}(elem.xq,elem.yq,AssembInfo.t);

        MxxSource = Data.sourceMxx_el{1}(elem.xq,elem.yq,AssembInfo.t);
        MxySource = Data.sourceMxy_el{1}(elem.xq,elem.yq,AssembInfo.t);
        MyxSource = Data.sourceMyx_el{1}(elem.xq,elem.yq,AssembInfo.t);
        MyySource = Data.sourceMyy_el{1}(elem.xq,elem.yq,AssembInfo.t);

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

