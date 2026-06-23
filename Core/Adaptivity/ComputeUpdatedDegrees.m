%> @file ComputeUpdatedDegrees.m
%> @author Mattia Corti, Caterina Leimer Saglio
%> @date 29 January 2026
%> @brief Compute updated degrees for the adaptivity
%>
%==========================================================================
%> @section classComputeUpdatedDegrees description
%==========================================================================
%> @brief Compute updated degrees for the adaptivity
%
%> @param Data        Struct with problem's data
%> @param mesh        Mesh struct (region+neighbor)
%> @param femregion   Finite Element struct (see CreateDOF.m)
%> @param Solution    Structure containing the numerical solution
%> @param AdIts       Iteration of p-adaptvity cycle
%>
%> @retval Data         Struct with updated data (computed the threshold)
%> @retval femregionNew Finite Element struct with updated degrees
%>                   
%==========================================================================

function [Data, femregionNew] = ComputeUpdatedDegrees(Data, mesh, femregion, Solution, AdIts)
        
    % Computation of centroids of the k-means clustering for tau threshold value
    if AdIts == 0
        [~,centroids] = kmeans(Solution.Indicator.tau,2);
        Data.tauThreshold = min(centroids);
    end

    % Update the degree according to the adaptation function
    degreeUpdated = Data.AdaptFunc(Solution.Indicator.tau/Data.tauThreshold);
    
    %% Create femregion with new degree distribution
    Data.degree = (femregion.degree>1).*(degreeUpdated<femregion.degree).*(femregion.degree-1) ...
                        + (femregion.degree<Data.maxDegree).*(degreeUpdated>femregion.degree).*(femregion.degree+1) ... 
                        + (degreeUpdated==femregion.degree).*femregion.degree;

    [femregionNew] = CreateDOF(Data, mesh.region);

    femregionNew.degree_old = femregion.degree;
    femregionNew.nbases_old = femregion.nbases;
end