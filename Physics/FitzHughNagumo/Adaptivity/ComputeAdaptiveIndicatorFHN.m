%> @file ComputeAdaptiveIndicatorFHN.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 11 May 2026
%> @brief Compute a-posteriori error indicator for adaptivity
%>
%==========================================================================
%> @section classComputeAdaptiveIndicatorFHN description
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

function [Indicator] = ComputeAdaptiveIndicatorFHN(Data, mesh, femregion, Solution, t)
    
    Funcs.Preallocation    = @IndicatorPreallocationFHN;
    Funcs.VolumeAssemblyST = @ResidualIndicatorAssemblyFHN;
    Funcs.FacesAssembly    = @FacesIndicatorsAssemblyFHN;
    Funcs.FinalMatrices    = @FinalIndicatorFHN;

    AssembInfo.quadrature              = "ST";
    AssembInfo.assemblyvolume          = true;
    AssembInfo.assemblyfaces           = true;
    AssembInfo.assemblyinternalfaces   = true;
    AssembInfo.assemblytrilinearforms  = false;
    AssembInfo.computegradients        = true;
    AssembInfo.computelaplacian        = true;
    AssembInfo.computefacegradients    = true;

    AssembInfo.u_h   = Solution.t.u_h;
    AssembInfo.w_h   = Solution.t.w_h;
    AssembInfo.u_old = Solution.t_old.u_h;
    AssembInfo.w_old = Solution.t_old.w_h;
    AssembInfo.t     = t;
    
    if isfield(femregion,'degree_old')
        AssembInfo.Matrices_adapt_old = Solution.t.Indicator.Adaptivity;
        AssembInfo.ass_vol_vec  = (femregion.degree~=femregion.degree_old);

        el_const = find(femregion.degree==femregion.degree_old);
        
        for kk = el_const'
            AssembInfo.ass_vol_vec(kk) = (femregion.degree(kk)>1) || ...
                (norm(Solution.t.u_h(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk)))-Solution.t_old.u_h(sum(femregion.nbases(1:kk-1))+1:sum(femregion.nbases(1:kk))))>1e-6);
        end
        AssembInfo.ass_face_vec = AssembInfo.ass_vol_vec;    
    else
        AssembInfo.ass_vol_vec  = ones(femregion.nel,1);
        AssembInfo.ass_face_vec = ones(femregion.nel,1);
    end

    fprintf("\n - Number of elements in which the indicator is updated: %i\n", sum(AssembInfo.ass_face_vec))

    [Indicator] = Assembly(Data, mesh.neighbor, femregion, AssembInfo, Funcs);
   
end
