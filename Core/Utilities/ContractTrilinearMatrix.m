%> @file  ContractTrilinearMatrix.m
%> @author Mattia Corti
%> @date 15 April 2026
%> @brief Compute the a matrix as the contraction of a third-order tensor M
%> with a vector x.
%>
%==========================================================================
%> @section classContractTrilinearMatrix Class description
%==========================================================================
%> @brief Contraction of a third-order tensor M with a vector x
%>
%> @param M                Third order tensor to be contracted
%> @param x                Vector \f$x\f$ to multiply
%> @param nel              Number of mesh elements
%> @param nbases           Number of bases associated to a single element
%>
%> @retval Mx              Matrix obtained by the contraction M*x
%>
%==========================================================================


function [Mx] = EvalTrilinearMatrix(M, x, nel, nbases)
    
    index = [0, cumsum(nbases')];

    for ie = 1:nel
        M{ie} = reshape(M{ie}(1:nbases(ie)^2,1:nbases(ie))*x(index(ie)+1:index(ie+1)),[nbases(ie), nbases(ie)]);
    end

    Mx = matlab.internal.math.blkdiag(M{:});

end
