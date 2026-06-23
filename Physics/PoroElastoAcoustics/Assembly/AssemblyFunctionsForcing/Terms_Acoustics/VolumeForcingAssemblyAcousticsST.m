%> @file   VolumeForcingAssemblyAcousticsST.m
%> @author Ilario Mazzieri, Mattia Corti
%> @date   30 May 2026
%> @brief Assembly of the vectors associated with volume integrals for
%> acoustics equation (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyAcousticsST Class description
%==========================================================================
%> @brief            Assembly of the vectors associated with volume
%> integrals for acoustics equation (subtriangulation implementation).
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

function [Forcing] = VolumeForcingAssemblyAcousticsST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters     
        rho_a    = Data.rho_a{id}(elem.xq,elem.yq);
        fSource1 = Data.source_phi{1}(elem.xq,elem.yq);

        % Vector assembling
        Forcing.F1_loc(1:nbases,1) = (elem.dx .* rho_a .* elem.phiq)' * fSource1;
        
end

