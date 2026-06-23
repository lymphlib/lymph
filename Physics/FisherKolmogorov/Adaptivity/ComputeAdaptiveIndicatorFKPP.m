%> @file ComputeAdaptiveIndicatorFKPP.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 14 April 2026
%> @brief Compute a-posteriori error indicator for adaptivity
%>
%==========================================================================
%> @section classComputeAdaptiveIndicatorFKPP description
%==========================================================================
%> @brief Compute a-posteriori error indicator for adaptivity
%
%> @param Data      Struct with problem's data
%> @param mesh      Mesh struct (region+neighbor)
%> @param femregion Finite Element struct (see CreateDOF.m)
%> @param Solution  Solution structure
%> @param t	        Current time
%>
%> @retval Indicator  A-posteriori error indicator 
%>                   
%==========================================================================

function [Indicator] = ComputeAdaptiveIndicatorFKPP(Data, mesh, femregion, Solution, t)
    
    Funcs.Preallocation    = @IndicatorPreallocationFKPP;
    Funcs.VolumeAssemblyST = @ResidualIndicatorAssemblyFKPP;
    Funcs.FacesAssembly    = @FacesIndicatorsAssemblyFKPP;
    Funcs.FinalMatrices    = @FinalIndicatorFKPP;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = true;
    AssembInfo.computefacegradients    = true;

    AssembInfo.c_h   = Solution.c_h;
    AssembInfo.c_old = Solution.c_old;
    AssembInfo.t     = t;
    
    if isfield(femregion,'degree_old')
        AssembInfo.Matrices_adapt_old = Solution.Indicator.Adaptivity;
        AssembInfo.ass_vol_vec  = (femregion.degree~=femregion.degree_old);

        el_const = find(femregion.degree==femregion.degree_old);
        
        for kk = el_const'
            AssembInfo.ass_vol_vec(kk) = (femregion.degree(kk)>1) || ...
                (norm(Solution.c_h(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk)))-Solution.c_old(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk))))>1e-6);
        end
        AssembInfo.ass_face_vec = AssembInfo.ass_vol_vec;    
    else
        AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
        AssembInfo.ass_face_vec = ones(femregion.nel,1);
    end

    fprintf("\n - Number of elements in which the indicator is updated: %i\n", sum(AssembInfo.ass_face_vec))

    [Indicator] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
