%> @file ComputeAdaptiveIndicatorLaplacian.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 12 May 2026
%> @brief Compute a-posteriori error indicator for adaptivity
%>
%==========================================================================
%> @section classComputeAdaptiveIndicatorLaplacian description
%==========================================================================
%> @brief Compute a-posteriori error indicator for adaptivity
%
%> @param Data        Struct with problem's data
%> @param mesh        Mesh struct (region+neighbor)
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Solution    Solution structure
%>
%> @retval Indicator  A-posteriori error indicator 
%>                   
%==========================================================================

function [Indicator] = ComputeAdaptiveIndicatorLaplacian(Data, mesh, femregion, Solution)
    
    Funcs.Preallocation    = @IndicatorPreallocationLaplacian;
    Funcs.VolumeAssemblyST = @ResidualIndicatorAssemblyLaplacian;
    Funcs.FacesAssembly    = @FacesIndicatorsAssemblyLaplacian;
    Funcs.FinalMatrices    = @FinalIndicatorLaplacian;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = true;
    AssembInfo.computefacegradients    = true;

    AssembInfo.uh = Solution.U;
    
    if isfield(femregion,'degree_old')
        AssembInfo.Matrices_adapt_old = Solution.Indicator.Adaptivity;
        AssembInfo.ass_vol_vec  = (femregion.degree~=femregion.degree_old);
        
        el_const = find(femregion.degree==femregion.degree_old);
        
        for kk = el_const'
            AssembInfo.ass_vol_vec(kk) = (norm(Solution.U(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk)))-Solution.Uold(sum(femregion.nbases_old(1:kk-1))+1:sum(femregion.nbases_old(1:kk))))>1e-6);
        end
        AssembInfo.ass_face_vec = AssembInfo.ass_vol_vec;
    else
        AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
        AssembInfo.ass_face_vec = ones(femregion.nel,1);
    end

    fprintf("\n - Number of elements in which the indicator is updated: %i\n", sum(AssembInfo.ass_face_vec))

    [Indicator] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end