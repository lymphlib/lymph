%> @file ComputeAdaptiveIndicatorHeat.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 6 February 2026
%> @brief Compute a-posteriori error indicator for adaptivity
%>
%==========================================================================
%> @section classComputeAdaptiveIndicatorHeat description
%==========================================================================
%> @brief Compute a-posteriori error indicator for adaptivity
%
%> @param Data        Struct with problem's data
%> @param mesh        Mesh struct (region+neighbor)
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Solution    Solution structure
%> @param t	      Current time
%>
%> @retval Indicator  A-posteriori error indicator 
%>                   
%==========================================================================

function [Indicator] = ComputeAdaptiveIndicatorHeat(Data, mesh, femregion, Solution, t)
    
    Funcs.Preallocation    = @IndicatorPreallocationHeat;
    Funcs.VolumeAssemblyST = @ResidualIndicatorAssemblyHeat;
    Funcs.FacesAssembly    = @FacesIndicatorsAssemblyHeat;
    Funcs.FinalMatrices    = @FinalIndicatorHeat;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = true;
    AssembInfo.computefacegradients    = true;

    AssembInfo.uh    = Solution.u_h;
    AssembInfo.u_old = Solution.u_old;
    AssembInfo.t     = t;
    
    if isfield(femregion,'degree_old')
        AssembInfo.Matrices_adapt_old = Solution.Indicator.Adaptivity;
        AssembInfo.ass_vol_vec  = (femregion.degree~=femregion.degree_old);
        
        el_const = find(femregion.degree==femregion.degree_old);
        
        for kk = el_const'
            AssembInfo.ass_vol_vec(kk) = (norm(Solution.u_h(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk)))-Solution.u_old(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk))))>1e-6);
        end
        AssembInfo.ass_face_vec = AssembInfo.ass_vol_vec;    
    else
        AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
        AssembInfo.ass_face_vec = ones(femregion.nel,1);
    end

    fprintf("\n - Number of elements in which the indicator is updated: %i\n", sum(AssembInfo.ass_face_vec))

    [Indicator] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
