%> @file   VolumeForcingAssemblyFHNST.m
%> @author Mattia Corti
%> @date   11 May 2026
%> @brief Assembly of the forcing terms with volume integrals for FHN
%> problem (subtriangulation implementation).
%>
%==========================================================================
%> @section classVolumeForcingAssemblyFHNST Class description
%==========================================================================
%> @brief            Assembly of the forcing terms associated with volume 
%> integrals for FHN problem (subtriangulation implementation).
%>
%> @param Data       Struct with problem's data
%> @param Forcing    Matrices struct (to be modified)
%> @param elem       Element struct containing the quadrature nodes and bases
%> @param ie         Element number
%> @param id         Element physical id 
%> @param nbases     Number of FE bases
%> @param AssembInfo Information associated to the specific assembly call
%>
%> @retval Forcing   Forcing struct
%>                   
%==========================================================================

function [Forcing] = VolumeForcingAssemblyFHNST(Data, Forcing, elem, ie, id, nbases, AssembInfo)
        
        % Evaluation of physical parameters
        k       = Data.k{id}(elem.xq,elem.yq,AssembInfo.t);
        a       = Data.a{id}(elem.xq,elem.yq,AssembInfo.t);
        epsilon = Data.epsilon{id}(elem.xq,elem.yq,AssembInfo.t);
        gamma   = Data.gamma{id}(elem.xq,elem.yq,AssembInfo.t);
        beta    = Data.beta{id}(elem.xq,elem.yq,AssembInfo.t);
        Chi     = Data.Chi{id}(elem.xq,elem.yq,AssembInfo.t);

        % solution reconstruction
        u_loc  = elem.phiq * AssembInfo.u(elem.index);
        w_loc  = elem.phiq * AssembInfo.w(elem.index);

        f      = Data.source_f{1}(elem.xq,elem.yq,AssembInfo.t);

        % Forcing term assembly
        Forcing.Iext_loc(1:nbases,1)  = (elem.dx .* elem.phiq)' * f ;
        Forcing.F_loc(1:nbases,1)     = (elem.dx .* elem.phiq)' * ( - Chi.*k.*u_loc.*(u_loc - 1).*(u_loc - a) - Chi.*w_loc );
        Forcing.G_loc(1:nbases,1)     = (elem.dx .* elem.phiq)' * ( epsilon.*(beta.*u_loc - gamma.*w_loc) );
        
end
                    
